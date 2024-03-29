from MA import *
from MSV import *
import random
from sv_util.os_aligners import *
from bokeh.plotting import figure, show
from sv_util.bokeh_style_helper import *
from sv_util.settings import *
from bokeh.plotting import ColumnDataSource
from bokeh.layouts import column, row, grid
from bokeh.models.tools import HoverTool
from svs_hidden_to_aligners.read_simulator import disfigure, load_sampled_from_file

def choice_adj_size(l, total_len):
    x = random.randrange(total_len)
    idx = 0
    while x >= l[idx][0]:
        x -= l[idx][0]
        idx += 1
    return l[idx]

def create_reads(pack, size, amount, func_get_seeds_and_read, do_disfigure=False):
    lumper = SeedLumping(ParameterSetManager())
    read_by_name = ReadByName()
    genome_section_by_name = ReadByName()
    points_by_name = {}
    contigs = [(x, y) for x, y in zip(pack.contigLengths(), pack.contigStarts()) if x > size]
    total_len = sum(x for x, _ in contigs)
    def read_and_seeds():
        contig_len, contig_start = choice_adj_size(contigs, total_len)
        start = random.randrange(contig_len - size)
        genome_section = pack.extract_from_to(start+contig_start, start+size+contig_start)
        return func_get_seeds_and_read(genome_section, start + contig_start), genome_section, start + contig_start
    i, d, m, il, dl, _ = load_sampled_from_file(sv_hidden_err_distrib)
    for idx in range(amount):
        read = NucSeq("N")
        while 'n' in str(read) or 'N' in str(read):
            (points, read), genome_section, sec_offset = read_and_seeds()
        if do_disfigure:
            s, up_shifts, right_shifts = disfigure(str(read), 1, d, m, i, il, dl)
            read = NucSeq(s)
            points_2 = []
            for p in points:
                p = list(p)
                for sp, sl in right_shifts:
                    if p[0] >= sp + sl:
                        p[0] -= sl
                        pass
                    elif p[0] >= sp and p[0] < sp + sl:
                        # this point is deleted
                        p = None
                        break
                if not p is None:
                    for sp, sl in up_shifts:
                        if p[0] >= sp:
                            p[0] += sl
                    points_2.append(p)
        else:
            points_2 = points
        read.name = "read" + str(idx)
        genome_section.name = "read" + str(idx)
        read_by_name.append(read)
        genome_section_by_name.append(genome_section)
        points_by_name[read.name] = (points_2, sec_offset)
    return points_by_name, read_by_name, genome_section_by_name


def create_scattered_read(pack, amount, num_pieces, size_pieces, do_disfigure=False):
    lumper = SeedLumping(ParameterSetManager())
    read_by_name = ReadByName()
    genome_section_by_name = ReadByName()
    points_by_name = {}
    contigs = [(x, y) for x, y in zip(pack.contigLengths(), pack.contigStarts()) if x > size_pieces]
    total_len = sum(x for x, _ in contigs)
    i, d, m, il, dl, _ = load_sampled_from_file(sv_hidden_err_distrib)
    def read_and_points():
        points = []
        read = ""
        for idx in range(0, num_pieces):
            contig_size, contig_start = choice_adj_size(contigs, total_len)
            start = random.randrange(contig_size - size_pieces)
            points.append((len(read), contig_start + start, True))
            points.append((len(read)+size_pieces, contig_start + start + size_pieces, True))
            read += str(pack.extract_from_to(contig_start+start,contig_start+start+size_pieces))
        return points, NucSeq(read)

    for idx in range(amount):
        read = NucSeq("N")
        while 'n' in str(read) or 'N' in str(read):
            points, read = read_and_points()
        if do_disfigure:
            s, up_shifts, right_shifts = disfigure(str(read), 1, d, m, i, il, dl)
            read = NucSeq(s)
            points_2 = []
            for p in points:
                p = list(p)
                for sp, sl in right_shifts:
                    if p[0] >= sp + sl:
                        p[0] -= sl
                    elif p[0] >= sp and p[0] < sp + sl:
                        # this point is deleted
                        p = None
                        break
                if not p is None:
                    for sp, sl in up_shifts:
                        if p[0] >= sp:
                            p[0] += sl
                    points_2.append(p)
        else:
            points_2 = points
        read.name = "read" + str(idx)
        read_by_name.append(read)
        genome_section = NucSeq()
        genome_section.name = "read" + str(idx)
        genome_section_by_name.append(genome_section)
        points_by_name[read.name] = (points_2, 0)
    return points_by_name, read_by_name, genome_section_by_name


def compare(params, data, points_by_name, reads, pack_pledge, fm_index_pledge,
            unlock_targets=None, render_one=True, original_seeds_list=None, ignore_genome_offset=False):
    compare_module = CompareSeedSets(params)
    lumper = SeedLumping(params)
    get_rectangles = SvJumpsFromSeeds(params, pack_pledge.get())
    collector = NucSeqSeedCollector(params)

    if render_one:
        for idx in range(params.get_num_threads()):
            read = reads[idx].get()
            while not read is None:
                #lumped_g_t = lumper.execute( ground_truth[0].get(), read, pack_pledge.get())
                lumped_data = lumper.execute( data[idx].get(), read, pack_pledge.get() )
                #if not original_seeds_list is None:
                #    printer = SeedPrinter(params)
                #    rectangles = get_rectangles.cpp_module.execute_helper(original_seeds_list[idx].get(),
                #                                                          pack_pledge.get(), read)
                #    printer.execute( lumped_data, rectangles )
                #else:
                if True:
                    printer = SeedPointPrinter(params)
                    printer.execute( lumped_data, points_by_name[read.name][0] )
                    exit()
                UnLock(params, unlock_targets[idx]).execute( Container() )
                read = reads[idx].get()
    else:
        res = VectorPledge()
        for idx in range(params.get_num_threads()):
            #lumped_g_t = promise_me(lumper, ground_truth[idx], reads[idx], pack_pledge)
            lumped_data = promise_me(lumper, data[idx], reads[idx], pack_pledge)
            empty = promise_me(collector, reads[idx], lumped_data)
            if unlock_targets is None:
                unlock = comp
            else:
                unlock = promise_me(UnLock(params, unlock_targets[idx]), empty)
            res.append(unlock)

        res.simultaneous_get(params.get_num_threads())

    def matches(seed, point, sec_offset):
        q,r,f = point
        if ignore_genome_offset:
            r -= sec_offset
        def nearby_start(max_diff=25):
            return abs(q-seed.start) <= max_diff and abs(r-seed.start_ref) <= max_diff
        def nearby_end(max_diff=25):
            if seed.on_forward_strand:
                return abs(q-(seed.start+seed.size)) <= max_diff and abs(r-(seed.start_ref+seed.size)) <= max_diff
            else:
                return abs(q-(seed.start+seed.size)) <= max_diff and abs(r-(seed.start_ref-seed.size)) <= max_diff
        return seed.on_forward_strand == f and (nearby_start() or nearby_end())

    all_found = {}
    for name, point_values in points_by_name.items():
        all_found[name] = [False]*len(point_values[0])
    for read, seeds in collector.cpp_module.collection:
        for seed in seeds:
            sec_offset = points_by_name[read.name][1]
            for idx, point in enumerate(points_by_name[read.name][0]):
                if matches(seed, point, sec_offset):
                    all_found[read.name][idx] = True

    s = 0
    for point_values in all_found.values():
        all_true = True
        for found in point_values:
            if not found:
                all_true = False
        if all_true:
            s += 1

    #print("hits:", s)

    return s / len(points_by_name)


def compare_seeds(params, reads_by_name, points_by_name, fm_index, pack, mems=True, reseeding=True,
                 render_one=False):
    #params.by_name("Number of Threads").set(1)
    #params.by_name("Use all Processor Cores").set(False)

    splitter = NucSeqSplitter(params)
    lock = Lock(params)
    reads_by_name_pledge = Pledge()
    reads_by_name_pledge.set(reads_by_name)
    pack_pledge = Pledge()
    pack_pledge.set(pack)
    fm_index_pledge = Pledge()
    fm_index_pledge.set(fm_index)
    min_len = MinLength(params, params.by_name("Minimal Seed Size SV").get() + 1)

    if mems:
        seeding_module = MinimizerSeeding(params)
        seed_lumping = SeedLumping(params)
        seed_extension = SeedExtender(params)
    else:
        seeding_module = BinarySeeding(params)
        extract_seeds = ExtractSeeds(params)

    reads_vec = ContainerVectorNucSeq()
    for name, read in reads_by_name:
        reads_vec.append(read)
    reads_vec_pledge = Pledge()
    reads_vec_pledge.set(reads_vec)

    data = []
    reads = []
    unlock_targets = []
    original_seeds_list = []
    read = promise_me(splitter, reads_vec_pledge)
    for _ in range(params.get_num_threads()):
        locked_read = promise_me(lock, read)
        unlock_targets.append(locked_read)
        reads.append(locked_read)
        if mems:
            minimizers = promise_me(seeding_module, fm_index_pledge, locked_read, pack_pledge)
            seeds = promise_me(seed_lumping, minimizers, locked_read, pack_pledge)
            if reseeding:
                soc_module = StripOfConsiderationSeeds(params)
                soc_filter = GetAllFeasibleSoCsAsSet(params)
                recursive_reseeding = RecursiveReseedingSoCs(params, pack)
                socs = promise_me(soc_module, seeds, locked_read, pack_pledge)
                soc_filtered_seeds = promise_me(soc_filter, socs)
                seeds = promise_me(recursive_reseeding, soc_filtered_seeds, pack_pledge, locked_read)
                original_seeds_list.append(seeds)
            else:
                seeds = promise_me(seed_extension, seeds, locked_read, pack_pledge)
                original_seeds_list = None
        else:
            segments = promise_me(seeding_module, fm_index_pledge, locked_read)
            seeds = promise_me(extract_seeds, segments, fm_index_pledge, locked_read)
            original_seeds_list.append(seeds)
            if reseeding:
                recursive_reseeding = RecursiveReseedingSegments(params, pack)
                seeds = promise_me(recursive_reseeding, segments, pack_pledge, fm_index_pledge, locked_read)
        data.append(seeds)

    return compare(params, data, points_by_name, reads, pack_pledge, fm_index_pledge,
                   unlock_targets, render_one, original_seeds_list=original_seeds_list)

def compare_nw(params, reads_by_name, genome_section_by_name, points_by_name, pack, render_one=False):
    #params.by_name("Number of Threads").set(1)
    #params.by_name("Use all Processor Cores").set(False)

    splitter = NucSeqSplitter(params)
    lock = Lock(params)
    genome_section_by_name_pledge = Pledge()
    genome_section_by_name_pledge.set(genome_section_by_name)
    pack_pledge = Pledge()
    pack_pledge.set(pack)

    reads_vec = ContainerVectorNucSeq()
    for name, read in reads_by_name:
        reads_vec.append(read)
    reads_vec_pledge = Pledge()
    reads_vec_pledge.set(reads_vec)

    nw_alignment = NWAlignment(params)
    alignment_to_seeds = AlignmentToSeeds(params)
    read_by_name = GetReadByReadName(params)

    data = []
    reads = []
    unlock_targets = []
    read = promise_me(splitter, reads_vec_pledge)
    for _ in range(params.get_num_threads()):
        locked_read = promise_me(lock, read)
        unlock_targets.append(locked_read)
        reads.append(locked_read)
        genome_section_pl = promise_me(read_by_name, locked_read, genome_section_by_name_pledge)
        alignment_pl = promise_me(nw_alignment, locked_read, genome_section_pl)
        seeds_pl = promise_me(alignment_to_seeds, alignment_pl, pack_pledge)
        data.append(seeds_pl)

    return compare(params, data, points_by_name, reads, pack_pledge, None, unlock_targets, render_one,
                    ignore_genome_offset=True)

def compare_alignment(params, reads_by_name_pledge, points_by_name, alignments, pack_pledge, fm_index_pledge,
                      unlock_targets=None, render_one=False):
    align_to_seeds = AlignmentToSeeds(params)
    get_read_by_name = GetReadByName(params)
    data = []
    reads = []
    for idx in range(params.get_num_threads()):
        alignment_seeds = promise_me(align_to_seeds, alignments[idx], pack_pledge)
        data.append(alignment_seeds)
        read = promise_me(get_read_by_name, alignments[idx], reads_by_name_pledge)
        reads.append(read)

    return compare(params, data, points_by_name, reads, pack_pledge, fm_index_pledge,
                   unlock_targets, render_one)



def compare_alignment_from_file_queue(params, reads_by_name, points_by_name, pack, fm_index, queue_pledge,
                                      render_one=False):
    file_reader = None # no op
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

    return compare_alignment(params, reads_by_name_pledge, points_by_name, alignments, pack_pledge, fm_index_pledge,
                             locked_files, render_one)


def compare_alignment_from_file_paths(params, reads_by_name, points_by_name, pack, fm_index, file_paths,
                                      render_one=False):
    if file_paths is None:
        return None
    file_queue = FileQueue()
    tot_size = 0
    for string in file_paths:
        if not os.path.exists(string):
            return float("NaN")
        tot_size += os.stat(string).st_size
        file_queue.add(FileStreamFromPath(string))
    if tot_size == 0:
        return 0
    queue_pledge = Pledge()
    queue_pledge.set(file_queue)
    return compare_alignment_from_file_queue(params, reads_by_name, points_by_name, pack, fm_index, queue_pledge,
                                             render_one)

def create_alignment(read_by_name, aligner, sam_name):
    reads_path = sv_hidden_to_aligners_data_dir + "/reads/" + sam_name + ".fasta"
    with open(reads_path, 'w') as fasta_file:
        for name, read in read_by_name:
            fasta_file.write(">" + name + "\n")
            fasta_file.write(str(read) + "\n")

    json_dict = {"reference_path":human_genome_dir}
    read_json = {"technology":"pb", "name":"n/a", "fasta_file":reads_path}
    path_sam = sv_hidden_to_aligners_data_dir + "/sam/" + sam_name + ".sam"
    aligner(read_json, path_sam, json_dict)

    return [path_sam]

def no_sv(genome_section, ref_start):
    seeds = Seeds()
    seeds.append(Seed(0, len(genome_section), ref_start, True))
    return seeds, genome_section

