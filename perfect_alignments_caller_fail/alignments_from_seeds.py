from MA import *
from MSV import *
import random
from sv_util.os_aligners import *
from bisect import bisect_right

def seeds_to_alignments(seeds, start, end, ref_pack):
    # Find rightmost seed less than or equal to start
    ret = []
    for seed in seeds:
        if seed.on_forward_strand:
            ret.append(Alignment(seed.start_ref, seed.start))
        else:
            ret.append(Alignment(ref_pack.pos_to_rev_strand(seed.start_ref), seed.start))
        ret[-1].append(MatchType.match, seed.size)
        ret[-1].mapping_quality = 60/254
    ret.sort(key=lambda x: x.get_score(), reverse=True)
    for alignment in ret[1:]:
        alignment.supplementary = True
    return ret


def crop_seeds(seeds, insertions, read_start, read_end):
    # Find rightmost seed less than or equal to start
    idx = bisect_right([seed.start for seed in seeds], read_start) - 1
    ret_seeds = Seeds()
    ret_ins = []
    q_pos = 0

    if idx < len(seeds) and seeds[idx].start + seeds[idx].size <= read_start:
        idx+=1
        s = read_end
        if idx < len(seeds) and seeds[idx].start < read_end:
            s = seeds[idx].start
        ins = str(insertions[idx-1])[-(s-read_start):]
        ret_seeds.append(Seed(q_pos, 0, s, seeds[idx].on_forward_strand))
        q_pos += len(ins)
        ret_ins.append(NucSeq(ins))
    while idx < len(seeds) and seeds[idx].start < read_end:
        s = max(seeds[idx].start, read_start)
        e = min(seeds[idx].start + seeds[idx].size, read_end)
        l = e-s
        if seeds[idx].on_forward_strand:
            s_r = seeds[idx].start_ref + s - seeds[idx].start
        else:
            s_r = seeds[idx].start_ref - s + seeds[idx].start
        if l > 0:
            ret_seeds.append(Seed(q_pos, l, s_r, seeds[idx].on_forward_strand))
            q_pos += l
        ret_ins.append(NucSeq(str(insertions[idx])[:read_end-e]))
        q_pos += min(len(insertions[idx]), read_end-e)
        idx += 1

    return ret_seeds, ret_ins


def rev_comp(read, seeds):
    read.complement()
    read.reverse()
    for seed in seeds:
        if seed.on_forward_strand:
            seed.start_ref += seed.size
        else:
            seed.start_ref -= seed.size
        seed.on_forward_strand = not seed.on_forward_strand
        seed.start = len(read) - (seed.start + seed.size)
    return read, seeds

def alignments_from_seeds(seeds, insertions, ref_pack, size, amount, paired_dist=None):
    ret = []
    y_end = seeds[-1].start + seeds[-1].size
    for idx in range(amount):
        def append_single():
            read_start = random.randrange(0, y_end - size)
            read_end = read_start + size
            c_seeds, c_ins = crop_seeds(seeds, insertions, read_start, read_end)
            read = reconstruct_sequenced_genome([("read_" + str(idx), c_seeds, c_ins)],
                                                           ref_pack).extract_forward_strand()
            read.name = "read_" + str(idx)
            read.id = idx
            if random.choice([True, False]):
                read, c_seeds = rev_comp(read, c_seeds)
            ret.append( (read, seeds_to_alignments(c_seeds, read_start, read_end, ref_pack)) )
        def append_paired():
            total_size = size * 2 + paired_dist
            read_start = random.randrange(0, y_end - total_size)
            read_end = read_start + total_size
            c_seeds, c_ins = crop_seeds(seeds, insertions, read_start, read_start + size)
            c_seeds_2, c_ins_2 = crop_seeds(seeds, insertions, read_end - size, read_end)
            read = reconstruct_sequenced_genome([("read_prim_" + str(idx*2), c_seeds, c_ins)],
                                                            ref_pack).extract_forward_strand()
            read.name = "read_prim_" + str(idx*2)
            read.id = idx*2
            read_2 = reconstruct_sequenced_genome([("read_mate_" + str(idx*2+1), c_seeds_2, c_ins_2)],
                                                              ref_pack).extract_forward_strand()
            read_2.name = "read_mate_" + str(idx*2+1)
            read_2.id = idx*2+1

            # check reads
            if len(read) != size:
                print(len(read), size)
                print(read)
                print(c_seeds)
                print(read_start, read_end)
                assert False
            if len(read_2) != size:
                print(len(read_2), size)
                print(read_2)
                print(c_seeds_2)
                print(read_start, read_end)
                assert False

            if random.choice([True, False]):
                read, c_seeds = rev_comp(read, c_seeds)
            else:
                read_2, c_seeds_2 = rev_comp(read_2, c_seeds_2)

            alignments = seeds_to_alignments(c_seeds, read_start, read_start + size, ref_pack)
            alignments_2 = seeds_to_alignments(c_seeds_2, read_end - size, read_end, ref_pack)
            for alignment in alignments:
                alignment.stats.first = True
                if len(alignments_2) > 0:
                    alignment.set_other(alignments_2[0]) # the first alignment in the list is the primary one
            for alignment in alignments_2:
                alignment.stats.first = False
                if len(alignments) > 0:
                    alignment.set_other(alignments[0]) # the first alignment in the list is the primary one
            alignments.extend(alignments_2)

            ret.append( ((read, read_2), alignments) )
        if paired_dist is None:
            append_single()
        else:
            append_paired()
    return ret

def alignment_to_file(alignments_list, sam_file_path, ref_pack, paired=False):
    params = ParameterSetManager()
    params.by_name("Emulate NGMLR's tag output").set(True)
    params.by_name("Soft clip").set(True)
    if paired:
        file_writer = PairedFileWriter(params, sam_file_path + ".sam", ref_pack)
        for (read_a, read_b), alignments in alignments_list:
            file_writer.execute(read_a, read_b, AlignmentVector(alignments), ref_pack)
    else:
        file_writer = FileWriter(params, sam_file_path + ".sam", ref_pack)
        for read, alignments in alignments_list:
            file_writer.execute(read, AlignmentVector(alignments), ref_pack)
    file_writer.cpp_module.close()
    sam_to_bam(sam_file_path)

def read_to_file(alignments_list, file_path):
    with open(file_path + ".fasta", "w") as fasta_file:
        for read, alignments in alignments_list:
            fasta_file.write(">" + str(read.name) + "\n")
            fasta_file.write(str(read) + "\n")
