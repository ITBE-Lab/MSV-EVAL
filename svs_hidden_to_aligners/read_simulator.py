import random
import json
import pysam
import numpy

##
# @brief apply modifications
# @details
# revcomp, deletions, mutations, insertions
# the order is important.
# code makes sure that the same nucleotide is never modified twice
# also indels must be at least one nuc apart from each other...
#
def disfigure(q, prob_modifier=0.1, del_prob=0.0277, mut_prob=0.0382, ins_prob=0.0651, ins_sizes=[1], 
              del_sizes=[1]):
    #deletions
    right_shifts = []
    pos = 0
    while pos < len(q):
        if random.random() < del_prob * prob_modifier:
            d_size = random.choice(del_sizes)
            right_shifts.append((pos, d_size))
            q = q[:pos] + q[pos + d_size:]
            pos += 1
        pos += 1

    # mutations
    pos = 0
    while pos < len(q):
        if random.random() < mut_prob * prob_modifier:
            q = q[:pos] + random.choice([c for c in "acgt" if c != q[pos].lower()]) + q[pos+1:]
        pos += 1

    # insertions
    pos = 0
    up_shifts = []
    while pos < len(q):
        if random.random() < ins_prob * prob_modifier:
            ins_size = random.choice(ins_sizes)
            up_shifts.append((pos,ins_size))
            ins_str = ""
            for _ in range(ins_size):
                ins_str += random.choice("acgt")
            q = q[:pos] + ins_str + q[pos:]
            pos += ins_size + 1
        else:
            pos += 1
    return q, up_shifts, right_shifts

def load_sampled_from_file(file):
    def _decode(o):
        if isinstance(o, str) or isinstance(o, unicode):
            try:
                return float(o)
            except ValueError:
                return o
        elif isinstance(o, dict):
            return {_decode(k): _decode(v) for k, v in o.items()}
        elif isinstance(o, list):
            return [_decode(v) for v in o]
        else:
            return o
    with open(file, "r") as f:
        json_file = json.loads(f.read(), object_hook=_decode)
        ins_prob, del_prob, mut_prob, ins_length_distrib, del_length_distrib, read_length_distrib = json_file
    return ins_prob, del_prob, mut_prob, ins_length_distrib, del_length_distrib, read_length_distrib

def get_cigars(sam_file):
    for alignment in pysam.AlignmentFile(sam_file):
        if not alignment.cigartuples is None:
            tag = None
            if alignment.has_tag("MD"):
                tag = alignment.get_tag("MD")
            yield (alignment.query_name, alignment.cigartuples, tag, alignment.query_length)

def get_characteristics_from_cigar(ins_length_distrib, del_length_distrib, cigar, tag):
    num_nuc = 0
    num_ins = 0
    num_del = 0
    num_ins_nuc = 0
    num_del_nuc = 0
    num_mm = 0
    bHaveM = False
    for operation, length in cigar:
            ## 0 -> M  BAM_CMATCH
            ## 1 -> I  BAM_CINS
            ## 2 -> D  BAM_CDEL
            ## 3 -> N  BAM_CREF_SKIP
            ## 4 -> S  BAM_CSOFT_CLIP
            ## 5 -> H  BAM_CHARD_CLIP
            ## 6 -> P  BAM_CPAD
            ## 7 -> =  BAM_CEQUAL
            ## 8 -> X  BAM_CDIFF
            ## 9 -> B  BAM_CBACK
            
            # match (=)
            if operation == 7: 
                num_nuc += length
            # insertion (I)
            elif operation == 1: 
                ins_length_distrib.append(length)
                num_ins += 1
                num_ins_nuc += length
            # deletion (D)
            elif operation == 2: 
                del_length_distrib.append(length)
                num_del += 1
                num_del_nuc += length
                num_nuc += length
            # mismatch (X)
            elif operation == 8: 
                num_mm += length
                num_nuc += length
            elif operation == 0:
                bHaveM = True
                num_nuc += length
            elif operation != 4 and operation != 5:
                print("got symbol:", operation, "in sam output")
    if bHaveM:
        if not tag is None and num_mm == 0:
            for index in range(len(tag)):
                if tag[index] in UPPERCASE_LETTERS:
                    if index == 0 or tag[index - 1] in NUMBERS or tag[index - 1] in UPPERCASE_LETTERS:
                        num_mm += 1
        else:
            print("cigar symbol M in alignment but", tag is None, num_mm == 0)
    return num_nuc, num_ins, num_del, num_mm, num_ins_nuc, num_del_nuc

def sample_distrib_from_fasta(
            sam_files = [""],
            out_file = "sampled_err_rate.json",
            num_reads = 1000000
        ):
    read_length_distrib = []
    print("Sampling read characteristics from sam...")
    num_nuc = 0
    num_ins = 0
    num_del = 0
    num_ins_nuc = 0
    num_del_nuc = 0
    num_mm = 0
    ins_length_distrib = []
    del_length_distrib = []
    cigar_list = []
    for sam_file in sam_files:
        print("\t" + sam_file)
        for idx, (_, cigar, tag, q_len) in enumerate(get_cigars(sam_file)):
            if idx > num_reads:
                break
            if idx % 10000 == 0:
                print(idx, "/", num_reads, "=", int(100*idx/num_reads), "%")
            read_length_distrib.append(q_len)
            num_nuc_, num_ins_, num_del_, num_mm_, num_ins_nuc_, num_del_nuc_ = get_characteristics_from_cigar(
                ins_length_distrib, del_length_distrib, cigar, tag)
            num_nuc += num_nuc_
            num_ins += num_ins_
            num_del += num_del_
            num_ins_nuc += num_ins_nuc_
            num_del_nuc += num_del_nuc_
            num_mm += num_mm_
    ins_prob = num_ins / num_nuc
    del_prob = num_del / num_nuc
    mut_prob = num_mm / num_nuc
    print("done")
    print("probabilities:")
    print("\tinsertion:", ins_prob)
    print("\tdeletion:", del_prob)
    print("\tmutation:", mut_prob)
    if mut_prob == 0:
        print("WARNING: aligner did not output any mutations...")
    print("indel sizes:")
    print(
            "\tinsertion [50% 98% 100%]:", 
            [numpy.percentile(ins_length_distrib, x) for x in [50, 98, 100]],
            "adjusted by length:",
            num_ins_nuc / num_nuc
        )
    print(
            "\tdeletion [50% 98% 100%]:", 
            [numpy.percentile(del_length_distrib, x) for x in [50, 98, 100]],
            "adjusted by length:",
            num_del_nuc / num_nuc
        )
    print(
            "read length [50% 98% 100%]:", 
            [numpy.percentile(read_length_distrib, x) for x in [50, 98, 100]]
        )
    with open(out_file, "w") as f:
        json.dump([
                ins_prob,
                del_prob,
                mut_prob,
                ins_length_distrib,
                del_length_distrib,
                read_length_distrib
            ], f)
    print("\tdone")
# from svs_hidden_to_aligners.read_simulator import *
# sample_distrib_from_fasta(["/MAdata/sv_caller_analysis/svs_hidden_to_aligners/public_revio_2022Q4_HG002-rep1/HG002.m84011_220902_175841_s1.GRCh38.bam"], "HG002.m84011_220902_175841_s1.GRCh38.sampled_errors.json")