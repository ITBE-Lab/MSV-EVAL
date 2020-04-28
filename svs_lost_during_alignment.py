from MA import *
from MSV import *
import random
from os_aligners import *
from bokeh.plotting import figure, show
from bokeh_style_helper import *
from bokeh.plotting import ColumnDataSource
from bokeh.layouts import column, row, grid
from bokeh.models.tools import HoverTool
from seed_printer import *

genome_dir = "/MAdata/genome/human/GRCh38.p12"
data_dir = "/MAdata/sv_lost_during_alignment"

def choice_adj_size(l, total_len):
    x = random.randrange(total_len)
    idx = 0
    while x >= l[idx][0]:
        x -= l[idx][0]
        idx += 1
    return l[idx]

def create_reads(pack, size, amount, func_get_seeds_and_read):
    lumper = SeedLumping(ParameterSetManager())
    read_by_name = ReadByName()
    seeds_by_name = SeedsByName()
    contigs = [(x, y) for x, y in zip(pack.contigLengths(), pack.contigStarts()) if x > size]
    total_len = sum(x for x, _ in contigs)
    def read_and_seeds():
        contig_len, contig_start = choice_adj_size(contigs, total_len)
        start = random.randrange(contig_len - size)
        genome_section = pack.extract_from_to(start+contig_start, start+size+contig_start)
        return func_get_seeds_and_read(genome_section, start + contig_start)
    for idx in range(amount):
        read = NucSeq("N")
        while 'n' in str(read) or 'N' in str(read):
            seeds, read = read_and_seeds()
        read.name = "read" + str(idx)
        read_by_name.append(read)
        lumped_seeds = lumper.execute(seeds, read, pack)
        seeds_by_name.append(lumped_seeds, read.name)
    return seeds_by_name, read_by_name


only_scattered_sv = True
no_scattered_sv = False
def create_scattered_read(pack, size, amount, num_pieces, size_pieces):
    read_by_name = ReadByName()
    seeds_by_name = SeedsByName()
    main_contigs = [(x, y) for x, y in zip(pack.contigLengths(), pack.contigStarts()) if x > size]
    main_total_len = sum(x for x, _ in main_contigs)
    contigs = [(x, y) for x, y in zip(pack.contigLengths(), pack.contigStarts()) if x > size_pieces]
    total_len = sum(x for x, _ in contigs)
    def read_and_seeds():
        seeds = Seeds()
        read = ""
        main_contig_len, main_cont_start = choice_adj_size(main_contigs, main_total_len)
        main_start = random.randrange(main_contig_len - size)
        size_one = size//(num_pieces+1)
        for idx in range(0, num_pieces+1):
            if not only_scattered_sv:
                seeds.append(Seed(len(read), size_one,
                                  main_cont_start + main_start+idx*size_one, True))
            read += str(pack.extract_from_to(main_cont_start+main_start+idx*size_one,
                                             main_cont_start+main_start+(idx+1)*size_one))
            if idx < num_pieces:
                contig_size, contig_start = choice_adj_size(contigs, total_len)
                start = random.randrange(contig_size - size_pieces)
                if not no_scattered_sv:
                    seeds.append(Seed(len(read), size_pieces, contig_start + start, True))
                read += str(pack.extract_from_to(contig_start+start,contig_start+start+size_pieces))
        return seeds, NucSeq(read)

    for idx in range(amount):
        read = NucSeq("N")
        while 'n' in str(read) or 'N' in str(read):
            seeds, read = read_and_seeds()
        read.name = "read" + str(idx)
        read_by_name.append(read)
        seeds_by_name.append(seeds, read.name)
    return seeds_by_name, read_by_name


def compare(params, ground_truth, data, seeds_by_name_pledge, reads, pack_pledge, fm_index_pledge,
            unlock_targets=None, render_one=False, segments_list=None):
    compare_module = CompareSeedSets(params)
    lumper = SeedLumping(params)
    get_seed_set_comp = GetSeedSetCompByName(params)
    get_rectangles = SvJumpsFromSeeds(params, pack_pledge.get())

    if render_one:
        for idx in range(params.get_num_threads()):
            read = reads[idx].get()
            while not read is None:
                #lumped_g_t = lumper.execute( ground_truth[0].get(), read, pack_pledge.get())
                lumped_data = lumper.execute( data[idx].get(), read, pack_pledge.get())
                printer = SeedPrinter(params)
                if not segments_list is None:
                    rectangles = get_rectangles.cpp_module.execute_helper(segments_list[idx].get(), pack_pledge.get(),
                                                               fm_index_pledge.get(), read)
                    printer.execute( lumped_data, ground_truth[idx].get(), rectangles )
                else:
                    printer.execute( lumped_data, ground_truth[idx].get() )
                UnLock(params, unlock_targets[idx]).execute( Container() )
                read = reads[idx].get()
    else:
        res = VectorPledge()
        for idx in range(params.get_num_threads()):
            #lumped_g_t = promise_me(lumper, ground_truth[idx], reads[idx], pack_pledge)
            lumped_data = promise_me(lumper, data[idx], reads[idx], pack_pledge)
            seed_set_comp = promise_me(get_seed_set_comp, reads[idx], seeds_by_name_pledge)
            comp = promise_me(compare_module, ground_truth[idx], lumped_data, seed_set_comp)
            if unlock_targets is None:
                unlock = comp
            else:
                unlock = promise_me(UnLock(params, unlock_targets[idx]), comp)
            res.append(unlock)

        res.simultaneous_get(params.get_num_threads())

    return seeds_by_name_pledge.get().mergeAll()


def compare_seeds(params, reads_by_name, seeds_by_name, fm_index, pack, reseeding=True, render_one=False):
    splitter = NucSeqSplitter(params)
    lock = Lock(params)
    reads_by_name_pledge = Pledge()
    reads_by_name_pledge.set(reads_by_name)
    pack_pledge = Pledge()
    pack_pledge.set(pack)
    fm_index_pledge = Pledge()
    fm_index_pledge.set(fm_index)
    seeds_by_name_pledge = Pledge()
    seeds_by_name_pledge.set(seeds_by_name)

    seeding_module = BinarySeeding(params)
    extract_seeds = ExtractSeeds(params)
    get_seeds_by_name = GetSeedsByReadName(params)

    reads_vec = ContainerVectorNucSeq()
    for name, read in reads_by_name:
        reads_vec.append(read)
    reads_vec_pledge = Pledge()
    reads_vec_pledge.set(reads_vec)

    ground_truth = []
    data = []
    reads = []
    unlock_targets = []
    segments_list = []
    read = promise_me(splitter, reads_vec_pledge)
    for _ in range(params.get_num_threads()):
        locked_read = promise_me(lock, read)
        unlock_targets.append(locked_read)
        reads.append(locked_read)
        segments = promise_me(seeding_module, fm_index_pledge, locked_read)
        if reseeding:
            recursive_reseeding = RecursiveReseeding(params, pack)
            segments_list.append(segments)
            seeds = promise_me(recursive_reseeding, segments, pack_pledge, fm_index_pledge, locked_read)
        else:
            rectangles = None
            seeds = promise_me(extract_seeds, segments, fm_index_pledge, locked_read)
        data.append(seeds)
        ground_truth_seeds = promise_me(get_seeds_by_name, locked_read, seeds_by_name_pledge)
        ground_truth.append(ground_truth_seeds)

    return compare(params, ground_truth, data, seeds_by_name_pledge, reads, pack_pledge, fm_index_pledge,
                   unlock_targets, render_one, segments_list=segments_list)

def compare_alignment(params, reads_by_name_pledge, seeds_by_name, alignments, pack_pledge, fm_index_pledge,
                      unlock_targets=None, render_one=False):
    align_to_seeds = AlignmentToSeeds(params)
    get_seeds_by_name = GetSeedsByName(params)
    get_read_by_name = GetReadByName(params)
    seeds_by_name_pledge = Pledge()
    seeds_by_name_pledge.set(seeds_by_name)
    ground_truth = []
    data = []
    reads = []
    for idx in range(params.get_num_threads()):
        alignment_seeds = promise_me(align_to_seeds, alignments[idx], pack_pledge)
        data.append(alignment_seeds)
        ground_truth_seeds = promise_me(get_seeds_by_name, alignments[idx], seeds_by_name_pledge)
        ground_truth.append(ground_truth_seeds)
        read = promise_me(get_read_by_name, alignments[idx], reads_by_name_pledge)
        reads.append(read)

    return compare(params, ground_truth, data, seeds_by_name_pledge, reads, pack_pledge, fm_index_pledge,
                   unlock_targets, render_one)



def compare_alignment_from_file_queue(params, reads_by_name, seeds_by_name, pack, fm_index, queue_pledge,
                                      render_one=False):
    file_reader = SamFileReader(params)
    queue_picker = FilePicker(params)
    queue_placer = FileAlignmentPlacer(params)
    lock = Lock(params)
    reads_by_name_pledge = Pledge()
    reads_by_name_pledge.set(reads_by_name)
    pack_pledge = Pledge()
    pack_pledge.set(pack)
    fm_index_pledge = Pledge()
    fm_index_pledge.set(fm_index)

    alignments = []
    locked_files = []
    for _ in range(params.get_num_threads()):
        picked_file = promise_me(queue_picker, queue_pledge)
        locked_file = promise_me(lock, picked_file)
        alignment_ = promise_me(file_reader, locked_file, pack_pledge, reads_by_name_pledge)
        alignment = promise_me(queue_placer, alignment_, locked_file, queue_pledge)
        locked_files.append(locked_file)
        alignments.append(alignment)

    return compare_alignment(params, reads_by_name_pledge, seeds_by_name, alignments, pack_pledge, fm_index_pledge,
                             locked_files, render_one)


def compare_alignment_from_file_paths(params, reads_by_name, seeds_by_name, pack, fm_index, file_paths,
                                      render_one=False):
    if file_paths is None:
        return None
    file_queue = FileQueue()
    for string in file_paths:
        file_queue.add(FileStreamFromPath(string))
    queue_pledge = Pledge()
    queue_pledge.set(file_queue)
    return compare_alignment_from_file_queue(params, reads_by_name, seeds_by_name, pack, fm_index, queue_pledge,
                                             render_one)

def create_alignment(read_by_name, aligner, sam_name):
    reads_path = data_dir + "/reads/" + sam_name + ".fasta"
    with open(reads_path, 'w') as fasta_file:
        for name, read in read_by_name:
            fasta_file.write(">" + name + "\n")
            fasta_file.write(str(read) + "\n")

    json_dict = {"reference_path":genome_dir}
    read_json = {"technology":"pb", "name":"n/a", "fasta_file":reads_path}
    path_sam = data_dir + "/mm_sam/" + sam_name + ".sam"
    aligner(read_json, path_sam, json_dict)

    return [path_sam]

def no_sv(genome_section, ref_start):
    seeds = Seeds()
    seeds.append(Seed(0, len(genome_section), ref_start, True))
    return seeds, genome_section

