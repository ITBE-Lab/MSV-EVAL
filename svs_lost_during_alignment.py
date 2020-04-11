from MA import *
import random
from os_aligners import *

genome_dir = "/MAdata/genome/human/GRCh38.p12-chr1"
data_dir = "/MAdata/sv_lost_during_alignment"

def choice_adj_size(l, total_len):
    x = random.randrange(total_len)
    idx = 0
    while x >= len(l[idx][0]):
        x -= len(l[idx][0])
        idx += 1
    return l[idx]

def from_to(contig, f, t):
    s = ""
    for x in range(f, f + t):
        s += contig[x]
    return NucSeq(s)

def create_reads(pack, size, amount, func_get_seeds_and_read):
    read_by_name = ReadByName()
    seeds_by_name = SeedsByName()
    contigs = [(x, y) for x, y in zip(pack.contigSeqs(), pack.contigStarts()) if len(x) > size]
    total_len = sum(len(x) for x, _ in contigs)
    def read_and_seeds():
        contig, contig_start = choice_adj_size(contigs, total_len)
        start = random.randrange(len(contig) - size)
        genome_section = from_to(contig, start, size)
        return func_get_seeds_and_read(genome_section, start + contig_start)
    for idx in range(amount):
        seeds, read = read_and_seeds()
        if 'n' in str(read) or 'N' in str(read):
            continue
        read.name = "read" + str(idx)
        read_by_name.append(read)
        seeds_by_name.append(seeds, read.name)
    return seeds_by_name, read_by_name


def compare(params, ground_truth, data, reads, pack_pledge, unlock_targets=None):
    collect = CollectSeedSetComps(params)
    compare_module = CompareSeedSets(params)
    lumper = SeedLumping(params)

    res = VectorPledge()
    for idx in range(params.get_num_threads()):
        lumped_g_t = promise_me(lumper, ground_truth[idx], reads[idx], pack_pledge)
        lumped_data = promise_me(lumper, data[idx], reads[idx], pack_pledge)
        comp = promise_me(compare_module, lumped_g_t, lumped_data)
        empty = promise_me(collect, comp)
        if unlock_targets is None:
            unlock = empty
        else:
            unlock = promise_me(UnLock(params, unlock_targets[idx]), empty)
        res.append(unlock)

    res.simultaneous_get(params.get_num_threads())
    return collect.cpp_module.collection


def compare_seeds(params, reads_by_name, seeds_by_name, fm_index, pack):
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
    read = promise_me(splitter, reads_vec_pledge)
    for _ in range(params.get_num_threads()):
        locked_read = promise_me(lock, read)
        unlock_targets.append(locked_read)
        reads.append(locked_read)
        segments = promise_me(seeding_module, fm_index_pledge, locked_read)
        seeds = promise_me(extract_seeds, segments, fm_index_pledge, locked_read)
        data.append(seeds)
        ground_truth_seeds = promise_me(get_seeds_by_name, locked_read, seeds_by_name_pledge)
        ground_truth.append(ground_truth_seeds)

    return compare(params, ground_truth, data, reads, pack_pledge, unlock_targets)

def compare_alignment(params, reads_by_name_pledge, seeds_by_name, alignments, pack_pledge, unlock_targets=None):
    align_to_seeds = AlignmentToSeeds(params)
    get_seeds_by_name = GetSeedsByName(params)
    get_read_by_name = GetReadByName(params)
    seeds_by_name_pledge = Pledge()
    seeds_by_name_pledge.set(seeds_by_name)
    ground_truth = []
    data = []
    reads = []
    for idx in range(params.get_num_threads()):
        alignment_seeds = promise_me(align_to_seeds, alignments[idx])
        data.append(alignment_seeds)
        ground_truth_seeds = promise_me(get_seeds_by_name, alignments[idx], seeds_by_name_pledge)
        ground_truth.append(ground_truth_seeds)
        read = promise_me(get_read_by_name, alignments[idx], reads_by_name_pledge)
        reads.append(read)

    return compare(params, ground_truth, data, reads, pack_pledge, unlock_targets)



def compare_alignment_from_file_queue(params, reads_by_name, seeds_by_name, pack, queue_pledge):
    file_reader = SamFileReader(params)
    queue_picker = FilePicker(params)
    queue_placer = FileAlignmentPlacer(params)
    lock = Lock(params)
    reads_by_name_pledge = Pledge()
    reads_by_name_pledge.set(reads_by_name)
    pack_pledge = Pledge()
    pack_pledge.set(pack)

    alignments = []
    locked_files = []
    for _ in range(params.get_num_threads()):
        picked_file = promise_me(queue_picker, queue_pledge)
        locked_file = promise_me(lock, picked_file)
        alignment_ = promise_me(file_reader, locked_file, pack_pledge, reads_by_name_pledge)
        alignment = promise_me(queue_placer, alignment_, locked_file, queue_pledge)
        locked_files.append(locked_file)
        alignments.append(alignment)

    return compare_alignment(params, reads_by_name_pledge, seeds_by_name, alignments, pack_pledge, locked_files)


def compare_alignment_from_file_paths(params, reads_by_name, seeds_by_name, pack, file_paths):
    if file_paths is None:
        return None
    file_queue = FileQueue()
    for string in file_paths:
        file_queue.add(FileStreamFromPath(string))
    queue_pledge = Pledge()
    queue_pledge.set(file_queue)
    return compare_alignment_from_file_queue(params, reads_by_name, seeds_by_name, pack, queue_pledge)

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

only_sv = True
no_sv = False
def translocation(offset, sv_size, gap_size, genome_section, ref_start):
    seeds = Seeds()
    total_size = sv_size*2 + gap_size + offset
    g = str(genome_section)
    # before translocation
    if not only_sv:
        seeds.append(Seed(0, offset, ref_start, True))
    read = g[:offset] 

    # section a
    if not no_sv:
        seeds.append(Seed(offset + sv_size + gap_size, sv_size, ref_start + offset, True))
    read += g[offset + sv_size + gap_size:total_size] 

    # within translocation
    if not only_sv:
        seeds.append(Seed(offset + sv_size, gap_size, ref_start + offset + sv_size, True))
    read += g[offset + sv_size:offset + sv_size + gap_size]

    # section b
    if not no_sv:
        seeds.append(Seed(offset, sv_size, ref_start + offset + sv_size + gap_size, True))
    read += g[offset:offset + sv_size]

    # after tranlocation
    if not only_sv:
        seeds.append(Seed(total_size, len(genome_section) - total_size, ref_start + total_size, True))
    read += g[total_size:]

    return seeds, NucSeq(read)

def main():
    params = ParameterSetManager()

    pack = Pack()
    pack.load(genome_dir + "/ma/genome")
    fm_index = FMIndex()
    fm_index.load(genome_dir + "/ma/genome")

    seeds_by_name, read_by_name = create_reads(pack, 1000, 100, lambda x,y: translocation(100, 20, 10, x,y))
    path_sam = create_alignment(read_by_name, mm2, "mm2")
    print("Minimap 2 alignment:")
    comp = compare_alignment_from_file_paths(params, read_by_name, seeds_by_name, pack, path_sam)
    print("overlapped:", 100 * comp.nt_overlap / comp.nt_ground_truth, "% (nt)",
          100 * comp.amount_overlap / comp.amount_ground_truth, "% (seeds)")
    print("SMEMs:")
    comp = compare_seeds(params, read_by_name, seeds_by_name, fm_index, pack)
    print("overlapped:", 100 * comp.nt_overlap / comp.nt_ground_truth, "% (nt)",
          100 * comp.amount_overlap / comp.amount_ground_truth, "% (seeds)")

main()