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
        if len(ret) > 1:
            ret[-1].supplementary = True
    return ret

def contigs_from_seeds(seeds, ref_pack):
    ret = []
    for seed in seeds:
        start_of_contig = ref_pack.start_of_sequence_id( ref_pack.seq_id_for_pos( seed.start_ref ) )
        if seed.start_ref > 0 and seed.start_ref == start_of_contig:
            ret.append(seed.start)

    return ret

def crop_seeds(seeds, insertions, read_start, read_end):
    # Find rightmost seed less than or equal to start
    idx = bisect_right([seed.start for seed in seeds], read_start) - 1
    ret_seeds = Seeds()
    ret_ins = []
    q_pos = 0
    
    while seeds[idx].start < read_end:
        s = max(seeds[idx].start, read_start)
        e = min(seeds[idx].start + seeds[idx].size, read_end)
        l = e-s
        if seeds[idx].on_forward_strand:
            s_r = (s - seeds[idx].start) + seeds[idx].start_ref
        else:
            # @todo this is not working....
            s_r = (s - seeds[idx].start) + (seeds[idx].start_ref - seeds[idx].size )
        if l > 0:
            ret_seeds.append(Seed(q_pos, l, s_r, seeds[idx].on_forward_strand))
            q_pos += l
        ret_ins.append(insertions[idx])
        q_pos += len(insertions[idx])
        idx += 1

    return ret_seeds, ret_ins


def alignments_from_db(call_table, ref_pack, caller_run, size, amount, start=0, end=None):
    if end is None:
        end = ref_pack.unpacked_size_single_strand
    ret = []

    seeds, insertions = call_table.calls_to_seeds(ref_pack, caller_run)
    contig_starts = contigs_from_seeds(seeds, ref_pack)
    seq_genome_size = max(seed.start + seed.size for seed in seeds)
    for _ in range(amount):
        while True:
            # @note no reverse strand reads for now
            read_start = random.randrange(start, end - size)
            read_end = read_start + size
            read_contig_id = bisect_right(contig_starts, read_start) - 1
            c_s = contig_starts[read_contig_id]
            if read_contig_id + 1 < len(contig_starts):
                c_l = contig_starts[read_contig_id + 1] - c_s
            else:
                c_l = seq_genome_size
            if read_start + size <= c_s + c_l:
                # break only if read is not bridging
                # @note this affects the coverage at chromosome endpoints...
                # shouldn't matter though - we will place our SVs with margins to the endpoints.
                break
        c_seeds, c_ins = crop_seeds(seeds, insertions, read_start, read_end)
        read = call_table.reconstruct_sequenced_genome_from_seeds(c_seeds, c_ins, ref_pack).extract_forward_strand()
        ret.append( (read, seeds_to_alignments(c_seeds, read_start, read_end, ref_pack)) )
    return ret

def alignment_to_file(alignments_list, sam_file_path, ref_pack):
    params = ParameterSetManager()
    params.by_name("Emulate NGMLR's tag output").set(True)
    file_writer = FileWriter(params, sam_file_path, ref_pack)
    for read, alignments in alignments_list:
        file_writer.execute(read, AlignmentVector(alignments), ref_pack)
    file_writer.close()
    sam_to_bam(sam_file_path)

