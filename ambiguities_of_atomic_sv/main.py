from MSV import *
from MA import *
import random
from ambiguities_of_atomic_sv.alignments_from_seeds import *
from bokeh.io import output_file
from ambiguities_of_atomic_sv.os_sv_callers import *
from ambiguities_of_atomic_sv.vcf_interpreters import *
import os
from sv_util.os_aligners import *
from ambiguities_of_atomic_sv.compute_errors import *
from sv_util.settings import *


def random_nuc_seq(l):
    ret = ""
    for _ in range(l):
        ret += random.choice(['A', 'C', 'G', 'T'])
    return ret

def four_nested_svs_calls(l, j):
    return ([
        Seed(0, j*l, 0, True),
        Seed((j+0)*l, l, (j+8)*l, False),
        Seed((j+1)*l, l, (j+3)*l, False),
        Seed((j+3)*l, l, (j+2)*l, False),
        Seed((j+4)*l, l, (j+4)*l, False),
        Seed((j+5)*l, l, (j+7)*l, False),
        Seed((j+6)*l, l, (j+5)*l, False),
        Seed((j+7)*l, l, (j+1)*l, False),
        Seed((j+8)*l, j*l, (j+8)*l, True),
    ], "four_nested_svs_calls", ["", "", random_nuc_seq(l), "", "", "", "", "", ""])

def inversion_in_inversion(l, j):
    return ([
        Seed(0, j*l, 0, True),
        Seed((j+0)*l, l, (j+3)*l, False),
        Seed((j+1)*l, l, (j+1)*l, True),
        Seed((j+2)*l, l, (j+1)*l, False),
        Seed((j+3)*l, j*l, (j+3)*l, True),
    ], "inversion_in_inversion", None)

def inversion_in_inversion_2(l, j):
    return ([
        Seed(0, j*l, 0, True),
        Seed((j+0)*l, l, (j+1)*l, True),
        Seed((j+1)*l, l, (j+1)*l, False),
        Seed((j+2)*l, j*l, (j+2)*l, True),
    ], "inversion_in_inversion_2", None)

def insertion_in_inversion(l, j):
    return ([
        Seed(0, j*l, 0, True),
        Seed((j+0)*l, l, (j+2)*l, False),
        Seed((j+2)*l, l, (j+1)*l, False),
        Seed((j+3)*l, j*l, (j+2)*l, True),
    ], "insertion_in_inversion", ["", random_nuc_seq(l), "", ""])

def inversion(l, j):
    return ([
        Seed(0, j*l, 0, True),
        Seed((j+0)*l, l, (j+1)*l, False),
        Seed((j+1)*l, j*l, (j+1)*l, True),
    ], "inversion", None)

def inversion_in_translocation(l, j):
    return ([
        Seed(0, j*l, 0, True),
        Seed((j+0)*l, l, (j+2)*l, True),
        Seed((j+1)*l, l, (j+2)*l, False),
        Seed((j+2)*l, l, (j+0)*l, True),
        Seed((j+3)*l, j*l, (j+3)*l, True),
    ], "inversion_in_translocation", None)

def inversion_in_translocation_separate_breakends(l, j):
    return ([
        Seed(0, j*l, 0, True),
        Seed((j+0)*l, l, (j+4)*l, True),
        Seed((j+1)*l, l, (j+1)*l, True),
        Seed((j+2)*l, l, (j+3)*l, False),
        Seed((j+3)*l, l, (j+3)*l, True),
        Seed((j+4)*l, l, (j+0)*l, True),
        Seed((j+5)*l, j*l, (j+5)*l, True),
    ], "inversion_in_translocation_separate_breakends", None)

def inversions_in_duplication(l, j):
    return ([
        Seed(0, j*l, 0, True),
        Seed((j+0)*l, l, (j+1)*l, False),
        Seed((j+1)*l, l, (j+1)*l, True),
        Seed((j+2)*l, l, (j+0)*l, True),
        Seed((j+3)*l, l, (j+2)*l, False),
        Seed((j+4)*l, j*l, (j+2)*l, True),
    ], "inversions_in_duplication", None)

def duplication_in_inversion(l, j):
    return ([
        Seed(0, j*l, 0, True),
        Seed((j+0)*l, l, (j+2)*l, False),
        Seed((j+1)*l, l, (j+2)*l, True),
        Seed((j+2)*l, l, (j+1)*l, False),
        Seed((j+3)*l, j*l, (j+2)*l, True),
    ], "duplication_in_inversion", None)

def inversion_overlapping_duplication(l, j):
    return ([
        Seed(0, (j+1)*l, 0, True),
        Seed((j+1)*l, 2*l, (j+2)*l, False),
        Seed((j+3)*l, (j+1)*l, (j+1)*l, True),
    ], "inversion_overlapping_duplication", None)

def translocation_in_duplication(l, j):
    return ([
        Seed(0, (j+1)*l, 0, True),
        Seed((j+1)*l, l, (j+2)*l, True),
        Seed((j+2)*l, l, (j+1)*l, True),
        Seed((j+3)*l, (j+4)*l, (j+0)*l, True),
    ], "translocation_in_duplication", None)

def duplication_of_inversion(l, j):
    return ([
        Seed(0, (j+1)*l, 0, True),
        Seed((j+1)*l, l, (j+3)*l, False),
        Seed((j+2)*l, l, (j+0)*l, True),
        Seed((j+3)*l, l, (j+3)*l, False),
        Seed((j+4)*l, l, (j+1)*l, True),
        Seed((j+5)*l, (j+0)*l, (j+3)*l, True),
    ], "duplication_of_inversion_with_translocation", None)

def overlapping_inversions_in_duplication(l, j):
    return ([
        Seed(0, (j+0)*l, 0, True),
        Seed((j+0)*l, l, (j+1)*l, False),
        Seed((j+1)*l, l, (j+1)*l, True),
        Seed((j+2)*l, 2*l, (j+2)*l, False),
        Seed((j+4)*l, (j+0)*l, (j+2)*l, True),
    ], "overlapping_inversions_in_duplication", None)

def inverted_duplication(l, j):
    return ([
        Seed(0, (j+1)*l, 0, True),
        Seed((j+1)*l, l, (j+1)*l, False),
        Seed((j+2)*l, (j+0)*l, (j+1)*l, True),
    ], "inverted_duplication", None)

def duplicated_inversion(l, j):
    return ([
        Seed(0, (j+0)*l, 0, True),
        Seed((j+0)*l, l, (j+1)*l, False),
        Seed((j+1)*l, l, (j+1)*l, False),
        Seed((j+2)*l, (j+0)*l, (j+1)*l, True),
    ], "duplicated_inversion", None)

def inversion_after_duplication(l, j):
    return ([
        Seed(0, (j+0)*l, 0, True),
        Seed((j+0)*l, l, (j+1)*l, False),
        Seed((j+1)*l, (j+1)*l, (j+0)*l, True),
    ], "inversion_after_duplication", None)

def triple_inversion(l, j):
    return ([
        Seed(0, (j+0)*l, 0, True),
        Seed((j+0)*l, l, (j+1)*l, True),
        Seed((j+1)*l, l, (j+0)*l, True),
        Seed((j+2)*l, (j+0)*l, (j+2)*l, True),
    ], "triple_inversion", None)

def overlapping_inversions(l, j):
    return ([
        Seed(0, (j+0)*l, 0, True),
        Seed((j+0)*l, l, (j+2)*l, False),
        Seed((j+1)*l, l, (j+3)*l, False),
        Seed((j+2)*l, l, (j+0)*l, True),
        Seed((j+3)*l, (j+0)*l, (j+3)*l, True),
    ], "overlapping_inversions", None)

def overlapping_inversions_2(l, j):
    return ([
        Seed(0, (j+0)*l, 0, True),
        Seed((j+0)*l, l, (j+2)*l, False),
        Seed((j+1)*l, l, (j+3)*l, False),
        Seed((j+2)*l, l, (j+1)*l, False),
        Seed((j+3)*l, (j+1)*l, (j+3)*l, True),
    ], "overlapping_inversions_2", None)

def negative_control(l, j):
    return ([
        Seed(0, (j+0)*l, 0, True),
    ], "negative_control", None)


def insertion_in_inversion_breakend(l, j):
    return ([
        Seed(0, j*l, 0, True),
        Seed((j+1)*l, l, (j+1)*l, False),
        Seed((j+2)*l, j*l, (j+1)*l, True),
    ], "insertion_in_inversion_breakend", [random_nuc_seq(l), "", ""])


def insertion_in_inversion_breakend_2(l, j):
    return ([
        Seed(0, j*l, 0, True),
        Seed((j+0)*l, l, (j+1)*l, False),
        Seed((j+2)*l, j*l, (j+1)*l, True),
    ], "insertion_in_inversion_breakend_2", ["", random_nuc_seq(l), ""])


def deletion_in_inversion_breakend(l, j):
    return ([
        Seed(0, j*l, 0, True),
        Seed((j+0)*l, l, (j+2)*l, False),
        Seed((j+1)*l, j*l, (j+2)*l, True),
    ], "deletion_in_inversion_breakend", ["", "", ""])


def deletion_in_inversion_breakend_2(l, j):
    return ([
        Seed(0, j*l, 0, True),
        Seed((j+0)*l, l, (j+1)*l, False),
        Seed((j+1)*l, j*l, (j+2)*l, True),
    ], "deletion_in_inversion_breakend_2", ["", "", ""])

def five_nested_inversions(l, j):
    return ([
        Seed(0, j*l, 0, True),                  # 0
        Seed((j+0)*l, l, (j+9)*l, False),       # 1
        Seed((j+1)*l, l, (j+1)*l, True),        # 2
        Seed((j+2)*l, l, (j+7)*l, False),       # 3
        Seed((j+3)*l, l, (j+3)*l, True),        # 4
        Seed((j+4)*l, l, (j+5)*l, False),       # 5
        Seed((j+5)*l, l, (j+5)*l, True),        # 6
        Seed((j+6)*l, l, (j+3)*l, False),       # 7
        Seed((j+7)*l, l, (j+7)*l, True),        # 8
        Seed((j+8)*l, l, (j+1)*l, False),       # 9
        Seed((j+9)*l, j*l, (j+9)*l, True),      # 10
    ], "five_nested_inversions", ["", "", "", "", "", "", "", "", "", "", ""])

def four_nested_inversions(l, j):
    return ([
        Seed(0, j*l, 0, True),                  # 0
        Seed((j+0)*l, l, (j+7)*l, False),       # 1
        Seed((j+1)*l, l, (j+1)*l, True),        # 2
        Seed((j+2)*l, l, (j+5)*l, False),       # 3
        Seed((j+3)*l, l, (j+3)*l, True),        # 4
        Seed((j+4)*l, l, (j+3)*l, False),       # 5
        Seed((j+5)*l, l, (j+5)*l, True),        # 6
        Seed((j+6)*l, l, (j+1)*l, False),       # 7
        Seed((j+7)*l, j*l, (j+7)*l, True),      # 8
    ], "four_nested_inversions", ["", "", "", "", "", "", "", "", ""])

def sequenced_genome_repetition(l, j):
    return ([
        Seed(0, j*l, 0, True),                  # 0
        Seed((j+0)*l, l, (j+1)*l, False),       # 1
        Seed((j+1)*l, 2*l, (j+1)*l, True),        # 2
        Seed((j+3)*l, 2*l, (j+0)*l, True),       # 3
        Seed((j+5)*l, l, (j+3)*l, False),        # 4
        Seed((j+6)*l, j*l, (j+3)*l, True),       # 5
    ], "sequenced_genome_repetition", ["", "", "", "", "", "",])

def deletion_and_duplication(l, j):
    return ([
        Seed(0, j*l, 0, True),                  # 0
        Seed((j+0)*l, 3*l, (j+1)*l, True),       # 1
        Seed((j+3)*l, l, (j+1)*l, True),        # 2
        Seed((j+4)*l, l, (j+3)*l, True),       # 3
        Seed((j+5)*l, j*l, (j+5)*l, True),        # 4
    ], "deletion_and_duplication", ["", "", "", "", "", "",])

def aBCDBDE(l, j):
    return ([
        Seed(0, (j+1)*l, 0, True),                  # 0
        Seed((j+1)*l, (j+4)*l, (j+2)*l, True),       # 1
    ], "ABCDBDE", ["", ""])

def aCABC(l, j):
    return ([
        Seed(0, (j+0)*l, 0, True),                  # 0
        Seed((j+0)*l, (j+4)*l, (j+1)*l, True),      # 1
    ], "ACABC", ["", ""])

def two_inverisons_in_duplication(l, j):
    return ([
        Seed(0, (j+1)*l, 0, True),                  # 0
        Seed((j+1)*l, l, (j+2)*l, False),           # 1
        Seed((j+2)*l, l, (j+2)*l, True),            # 2
        Seed((j+3)*l, l, (j+0)*l, True),            # 3
        Seed((j+4)*l, (j+1)*l, (j+2)*l, True),      # 4
    ], "two_inverisons_in_duplication", ["", "", "", "", ""])

def duplicated_inversion_flanked(l, j):
    return ([
        Seed(0, (j+1)*l, 0, True),                  # 0
        Seed((j+1)*l, l, (j+2)*l, True),            # 1
        Seed((j+2)*l, l, (j+3)*l, False),           # 2
        Seed((j+3)*l, (j+1)*l, (j+4)*l, True),      # 3
    ], "duplicated_inversion_flanked", ["", "", "", ""])

def translocations_in_duplication(l, j):
    return ([
        Seed(0, (j+1)*l, 0, True),                  # 0
        Seed((j+1)*l, l, (j+3)*l, True),            # 1
        Seed((j+2)*l, l, (j+0)*l, True),            # 2
        Seed((j+3)*l, l, (j+2)*l, True),            # 3
        Seed((j+4)*l, l, (j+1)*l, True),            # 4
        Seed((j+5)*l, (j+1)*l, (j+3)*l, True),      # 5
    ], "translocations_in_duplication", ["", "", "", "", "", ""])

def del_dup_inv(l, j):
    return ([
        Seed(0, (j+1)*l, 0, True),                  # 0
        Seed((j+1)*l, l, (j+2)*l, True),            # 1
        Seed((j+2)*l, l, (j+3)*l, False),           # 2
        Seed((j+3)*l, (j+1)*l, (j+3)*l, True),      # 3
    ], "del_dup_inv", ["", "", "", ""])

def double_insertion_same(l, j):
    ins = random_nuc_seq(l)
    return ([
        Seed(0, (j+1)*l, 0, True),                  # 0
        Seed((j+2)*l, l, (j+2)*l, False),            # 1
        Seed((j+4)*l, (j+1)*l, (j+2)*l, True),      # 2
    ], "double_insertion_same", [ins, ins, ""])

def double_insertion_diff(l, j):
    return ([
        Seed(0, (j+1)*l, 0, True),                  # 0
        Seed((j+2)*l, l, (j+2)*l, False),            # 1
        Seed((j+4)*l, (j+1)*l, (j+2)*l, True),      # 2
    ], "double_insertion_diff", [random_nuc_seq(l), random_nuc_seq(l), ""])

def insertion_and_deletion(l, j):
    return ([
        Seed(0, (j+1)*l, 0, True),                  # 0
        Seed((j+2)*l, (j+1)*l, (j+2)*l, True),      # 1
    ], "insertion_and_deletion", [random_nuc_seq(l), ""])

def insertion_and_inversion(l, j):
    return ([
        Seed(0, (j+1)*l, 0, True),                  # 0
        Seed((j+2)*l, (j+1)*l, (j+2)*l, False),      # 1
    ], "insertion_and_inversion", [random_nuc_seq(l), ""])

def insertion(l, j):
    return ([
        Seed(0, (j+1)*l, 0, True),                  # 0
        Seed((j+2)*l, (j+1)*l, (j+1)*l, True),      # 1
    ], "insertion", [random_nuc_seq(l), ""])

def check_insertion(seq, ref, insertions, j, l):
    assert seq[:(j+1)*l] == ref[:(j+1)*l]
    assert seq[(j+1)*l:(j+2)*l] == insertions[0]
    assert seq[(j+2)*l:] == ref[(j+1)*l:]
    assert seq == ref[:(j+1)*l] + insertions[0] + ref[(j+1)*l:]

def duplication(l, j):
    return ([
        Seed(0, (j+1)*l, 0, True),                  # 0
        Seed((j+1)*l, (j+1)*l, (j+0)*l, True),      # 1
    ], "duplication", ["", ""])

def duplication_and_insertion(l, j):
    return ([
        Seed(0, (j+1)*l, 0, True),                  # 0
        Seed((j+2)*l, (j+1)*l, (j+0)*l, True),      # 1
    ], "duplication_and_insertion", [random_nuc_seq(l), ""])

def inversion_after_inversion(l, j):
    return ([
        Seed(0, (j+1)*l, 0, True),                  # 0
        Seed((j+1)*l, l, (j+2)*l, False),          # 1
        Seed((j+2)*l, l, (j+2)*l, True),           # 2
        Seed((j+3)*l, l, (j+4)*l, False),          # 3
        Seed((j+4)*l, (j+1)*l, (j+4)*l, True),      # 4
    ], "inversion_after_inversion", ["", "", "", "", ""])

#l = 1000
l = 100
j = 10
coverage = 100
read_size = 2*l-1
callers = [
    (sniffles, "sniffles", sniffles_interpreter, "single"),
    #(manta, "manta", manta_interpreter, "paired"),
    (delly, "delly", delly_interpreter, "paired"),
    (gridss, "gridss", gridss_interpreter, "paired"),
    #(pbSv, "pbSv", pb_sv_interpreter),
]
paired_dist = 100
paired_size = 250

def run_msv(pack, seeds, insertions):
    s = libMA.containers.Seeds(seeds)
    n = [NucSeq(x) for x in insertions]
    param = ParameterSetManager()
    param.by_name("Do Dummy Jumps").set(False)
    reconstructed_query_genome = reconstruct_sequenced_genome([("sequenced_genome", s, n)],
                                                              pack).extract_forward_strand()
    return SvJumpsFromExtractedSeeds(param, pack).execute(s, pack,
                                                          reconstructed_query_genome)


def get_ref_section(section_size, j, l):
    ref_section = "N"
    while 'n' in ref_section or 'N' in ref_section:
        offset = random.randrange(1000, chr1_len - section_size - 1000)
        ref_section = str(reference.extract_from_to(offset, offset+section_size))
    return ref_section, ""

def get_ABCDBDE_section(section_size, j, l):
    section_size = (2*j+5)*l
    ref_section = "N"
    while 'n' in ref_section or 'N' in ref_section:
        offset = random.randrange(1000, chr1_len - section_size - 1000)
        ref_section = str(reference.extract_from_to(offset, offset+section_size))
    a = ref_section[:j*l]
    b = ref_section[j*l:(j+1)*l]
    c = ref_section[(j+1)*l:(j+2)*l]
    d = ref_section[(j+2)*l:(j+3)*l]
    e = ref_section[(j+3)*l:(2*j+4)*l]
    return a + b + c + d + b + d + e, ""

def get_ABCABC_section(section_size, j, l):
    section_size = (2*j+3)*l
    ref_section = "N"
    while 'n' in ref_section or 'N' in ref_section:
        offset = random.randrange(1000, chr1_len - section_size - 1000)
        ref_section = str(reference.extract_from_to(offset, offset+section_size))
    #print("get_ABCABC_section", j*l, (j+3)*l)
    f = ref_section[:j*l]
    abc = ref_section[j*l:(j+3)*l]
    e = ref_section[(j+3)*l:]
    return f + abc + abc + e, ""


if __name__ == "__main__":
    print("HAVE YOU ACTIVATED THE CONDA ENVIRONMENT?")
    reference = Pack()
    reference.load(human_genome_dir + "/ma/genome")

    chr1_len = reference.contigLengths()[0]

    svs = [
        #(negative_control, get_ref_section, None),
        #(inversion, get_ref_section, None),
        #(four_nested_svs_calls, get_ref_section, None),
        #(inversion_in_inversion, get_ref_section, None),
        (inversion_in_inversion_2, get_ref_section, None), # selected
        #(insertion_in_inversion, get_ref_section, None),
        #(inversion_in_translocation, get_ref_section, None),
        #(inversion_in_translocation_separate_breakends, get_ref_section, None),
        #(inversions_in_duplication, get_ref_section, None),
        #(duplication_in_inversion, get_ref_section, None),
        #(inversion_overlapping_duplication, get_ref_section, None), # selected
        #(translocation_in_duplication, get_ref_section, None),
        #(duplication_of_inversion, get_ref_section, None),
        #(overlapping_inversions_in_duplication, get_ref_section, None), # selected
        (inverted_duplication, get_ref_section, None), # selected
        (duplicated_inversion, get_ref_section, None), # selected
        (inversion_after_duplication, get_ref_section, None), # selected
        #(triple_inversion, get_ref_section, None),
        #(overlapping_inversions, get_ref_section, None),
        #(overlapping_inversions_2, get_ref_section, None),
        #(insertion_in_inversion_breakend, get_ref_section, None),
        #(insertion_in_inversion_breakend_2, get_ref_section, None),
        #(deletion_in_inversion_breakend, get_ref_section, None),
        #(deletion_in_inversion_breakend_2, get_ref_section, None),
        #(five_nested_inversions, get_ref_section, None),
        #(four_nested_inversions, get_ref_section, None),
        #(sequenced_genome_repetition, get_ref_section, None),
        #(deletion_and_duplication, get_ref_section, None),
        #(aBCDBDE, get_ABCDBDE_section, None),
        #(aCABC, get_ABCABC_section, None),
        #(two_inverisons_in_duplication, get_ref_section, None),
        #(duplicated_inversion_flanked, get_ref_section, None),
        #(translocations_in_duplication, get_ref_section, None),
        #(del_dup_inv, get_ref_section, None),
        #(double_insertion_same, get_ref_section, None),
        #(double_insertion_diff, get_ref_section, None),
        #(insertion_and_deletion, get_ref_section, None),
        #(insertion_and_inversion, get_ref_section, None),
        #(insertion, get_ref_section, check_insertion),
        #(duplication, get_ref_section, None),
        #(duplication_and_insertion, get_ref_section, None),
        (inversion_after_inversion, get_ref_section, None),
    ]

    output_file(ambiguities_of_atomic_sv_data_dir + "/bokeh_out_perfect_alignments_experiment.html")
    sets = []
    for sv_func, get_ref, check_genomes in svs:
        ref_section = "N"
        seeds, name, insertions = sv_func(l, j)
        section_size = seeds[-1].start_ref + seeds[-1].size
        ref_section, name_suffix = get_ref(section_size, j, l)
        name = name + name_suffix
        file_name = "perfect_alignments_experiment-" + name
        if insertions is None:
            insertions = [""] * len(seeds)

        reference_section = Pack()
        reference_section.append("chr1", "section of chr1", NucSeq(ref_section))

        if not check_genomes is None:
            [(seq_genome, _)] = alignments_from_seeds(seeds, insertions, reference_section,
                                                      seeds[-1].start + seeds[-1].size, 1, do_revcomp=False)
            check_genomes(str(seq_genome), str(ref_section), insertions, j, l)

        # create fake genome region
        with open(fasta_folder + file_name + ".genome.fasta", "w") as genome_out:
            genome_out.write(">chr1\n")
            genome_out.write(ref_section)
            genome_out.write("\n")
        if os.path.exists(fasta_folder + file_name + ".genome.fasta.fai"):
            os.system("rm " + fasta_folder + file_name + ".genome.fasta.*")
        os.system(sam_tools_pref + "faidx " + fasta_folder + file_name + ".genome.fasta")

        #for seed in seeds:
        #    seed.start_ref += offset

        num_reads = (coverage * section_size) // read_size
        # alignments with SVs
        alignments_list = alignments_from_seeds(seeds, insertions, reference_section, read_size, num_reads)

        num_reads_paired = (coverage * section_size) // (paired_size*2)
        # alignments with SVs
        alignments_list_paired = alignments_from_seeds(seeds, insertions, reference_section, paired_size,
                                                        num_reads_paired, paired_dist=paired_dist, do_revcomp=False)

        if False:
            seed_printer = SeedPrinter(ParameterSetManager(), "alignment seed", x_range=(offset, offset+section_size),
                                    y_range=(0, read_size), do_print=False)
            for _, alignments in alignments_list:
                alignment_seeds = Seeds()
                for alignment in alignments:
                    print(alignment.cigarString(reference_section, paired_size, False))
                    alignment_seeds.extend( alignment.to_seeds(reference_section) )
                seed_printer.execute(alignment_seeds)
                exit()


        sam_file_name = sam_folder + file_name
        sam_file_name_paired = sam_folder + file_name + "-paired"
        if False: # write alignments myself
            alignment_to_file(alignments_list[:1], sam_folder + "us", reference_section)

        alignment_to_file(alignments_list, sam_file_name, reference_section)
        alignment_to_file(alignments_list_paired, sam_file_name_paired, reference_section, paired=True)

        msv_entries = run_msv(reference_section, seeds, insertions)
        # other callers
        from_to_calls_lists = []
        with open(ambiguities_of_atomic_sv_data_dir + "/vcf_errors.log", "w") as error_file:
            print("dataset:", name)
            for caller, caller_name, interpreter, read_type in callers:
                vcf_file_path = vcf_folder + file_name + "-" + caller_name + ".vcf"
                if os.path.exists(vcf_file_path):
                    os.remove(vcf_file_path)
                caller( (sam_file_name if read_type=="single" else sam_file_name_paired) + ".sorted.bam",
                        vcf_file_path, fasta_folder + file_name + ".genome.fasta")
                if not os.path.exists( vcf_file_path ):
                    print("caller did not create calls:", caller_name)
                else:
                    print("caller did create calls:", caller_name)
                    from_to_calls_list = []
                    for call in vcf_parser(vcf_file_path):
                        a = interpreter(call, reference_section, error_file)
                        if a is None:
                            continue
                        for x in a:
                            from_to_calls_list.append(x)
                    from_to_calls_lists.append((caller_name, [*from_to_calls_list]))
        sets.append( (seeds, from_to_calls_lists, name, msv_entries) )
    render(sets, l)