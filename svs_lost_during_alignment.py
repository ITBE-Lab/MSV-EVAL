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


def compare(params, ground_truth, data, unlock_targets=None):
    collect = CollectSeedSetComps(params)
    compare_module = CompareSeedSets(params)

    res = VectorPledge()
    for idx in range(params.get_num_threads()):
        comp = promise_me(compare_module, ground_truth[idx], data[idx])
        empty = promise_me(collect, comp)
        if unlock_targets is None:
            unlock = empty
        else:
            unlock = promise_me(UnLock(params, unlock_targets[idx]), empty)
        res.append(unlock)

    print("simultaneous_get...")
    res.simultaneous_get(params.get_num_threads())
    print("simultaneous_get done")
    return collect.cpp_module.collection

def compare_alignment(params, reads_by_name, seeds_by_name, alignments, unlock_targets=None):
    ground_truth = []
    data = []
    align_to_seeds = AlignmentToSeeds(params)
    get_seeds_by_name = GetSeedsByName(params)
    seeds_by_name_pledge = Pledge()
    seeds_by_name_pledge.set(seeds_by_name)
    for idx in range(params.get_num_threads()):
        alignment_seeds = promise_me(align_to_seeds, alignments[idx])
        data.append(alignment_seeds)
        ground_truth_seeds = promise_me(get_seeds_by_name, alignments[idx], seeds_by_name_pledge)
        ground_truth.append(ground_truth_seeds)

    return compare(params, ground_truth, data, unlock_targets)

def compare_alignment_from_file_queue(params, reads_by_name, seeds_by_name, pack, queue_pledge):
    alignments = []
    locked_files = []
    file_reader = SamFileReader(params)
    queue_picker = FilePicker(params)
    queue_placer = FileAlignmentPlacer(params)
    lock = Lock(params)
    reads_by_name_pledge = Pledge()
    reads_by_name_pledge.set(reads_by_name)
    pack_pledge = Pledge()
    pack_pledge.set(pack)
    for _ in range(params.get_num_threads()):
        picked_file = promise_me(queue_picker, queue_pledge)
        locked_file = promise_me(lock, picked_file)
        alignment_ = promise_me(file_reader, locked_file, pack_pledge, reads_by_name_pledge)
        alignment = promise_me(queue_placer, alignment_, locked_file, queue_pledge)
        locked_files.append(locked_file)
        alignments.append(alignment)

    return compare_alignment(params, reads_by_name, seeds_by_name, alignments, locked_files)


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

def translocation(offset, sv_size, gap_size, genome_section, ref_start):
    seeds = Seeds()
    seeds.append(Seed(0, offset, ref_start, True))
    seeds.append(Seed(offset + sv_size, offset, ref_start + offset + sv_size, True))
    seeds.append(Seed(offset + sv_size + gap_size, offset, ref_start + offset + sv_size, True))
    seeds.append(Seed(offset + sv_size, offset, ref_start + offset + sv_size + gap_size, True))
    total_size = sv_size + gap_size * 2
    seeds.append(Seed(offset + total_size, len(genome_section) - total_size, ref_start + total_size, True))
    g = str(genome_section)
    read = g[:offset] 
    read += g[offset + sv_size + gap_size:total_size] 
    read += g[offset + sv_size:offset + sv_size + gap_size]
    read += g[offset:offset + sv_size]
    read += g[total_size:]
    return seeds, NucSeq(read)

def main():
    params = ParameterSetManager()

    pack = Pack()
    pack.load(genome_dir + "/ma/genome")

    print("creating reads...")
    seeds_by_name, read_by_name = create_reads(pack, 1000, 100, lambda x,y: translocation(100, 20, 10, x,y))
    print("creating alignments...")
    path_sam = create_alignment(read_by_name, mm2, "mm2")
    print("compating alignments...")
    comp = compare_alignment_from_file_paths(params, read_by_name, seeds_by_name, pack, path_sam)
    print("overlapped:", 100 * comp.nt_overlap / comp.nt_ground_truth, "% (nt)",
          100 * comp.amount_overlap / comp.amount_ground_truth, "% (seeds)")
    print("entropy:",100 * comp.nt_overlap / comp.nt_data, "% (nt)",
          100 * comp.amount_overlap / comp.amount_data, "% (seeds)")

main()