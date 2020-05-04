from MSV import *
from MA import *
import random
from perfect_alignments_caller_fail.alignments_from_db import *
from bokeh.io.output import reset_output
from caller_analysis.os_sv_callers import *
from caller_analysis.vcf_interpreters import *
import os
from sv_util.os_aligners import *

global_prefix = "/MAdata/sv_lost_during_calling/"
sam_folder = global_prefix + "sam/"
fasta_folder = global_prefix + "fasta/"
vcf_folder = global_prefix + "vcf/"
genome_dir = "/MAdata/genome/human"
db_prefix = "/MAdata/sv_datasets/"

def random_nuc_seq(l):
    ret = ""
    for _ in range(l):
        ret += random.choice(['a', 'c', 'g', 't'])
    return ret

#
# the genome reconstruction test ppt file contains the calls inserted here
#
def four_nested_svs_calls(db_conn, dataset_name, l, offset):
    JumpRunTable(db_conn)
    SvCallerRunTable(db_conn)
    get_inserter = GetCallInserter(ParameterSetManager(), db_conn, "four_nested_svs_calls",
                                   "the sv's that were simulated", -1)
    pool = PoolContainer(1, dataset_name)
    sv_inserter = get_inserter.execute(pool)

    sv_inserter.insert(SvCall(offset + l - 1, offset + 9*l - 1, 0, 0, True, 1000)) # a
    sv_inserter.insert(SvCall(offset + 4*l - 1, offset + 8*l, 0, 0, False, 1000)) # b

    insertion = SvCall(offset + 3*l - 1, offset + 3*l, 0, 0, False, 1000)
    insertion.inserted_sequence = NucSeq(random_nuc_seq(l))
    sv_inserter.insert(insertion) # c

    sv_inserter.insert(SvCall(offset + 5*l - 1, offset + 2*l, 0, 0, False, 1000)) # d
    sv_inserter.insert(SvCall(offset + 8*l - 1, offset + 4*l, 0, 0, False, 1000)) # e
    sv_inserter.insert(SvCall(offset + 6*l - 1, offset + 7*l, 0, 0, False, 1000)) # f
    sv_inserter.insert(SvCall(offset + 2*l - 1, offset + 5*l, 0, 0, False, 1000)) # g
    sv_inserter.insert(SvCall(offset + 9*l, offset + l, 0, 0, True, 1000)) # h

    sv_inserter.close(pool)
    return get_inserter.cpp_module.id

def inversion_in_inversion(db_conn, dataset_name, l, offset):
    JumpRunTable(db_conn)
    SvCallerRunTable(db_conn)
    get_inserter = GetCallInserter(ParameterSetManager(), db_conn, "inversion_in_inversion",
                                   "the sv's that were simulated", -1)
    pool = PoolContainer(1, dataset_name)
    sv_inserter = get_inserter.execute(pool)

    sv_inserter.insert(SvCall(offset + l - 1, offset + 4*l - 1, 0, 0, True, 1000)) # a
    sv_inserter.insert(SvCall(offset + 2*l, offset + 3*l, 0, 0, True, 1000)) # b
    sv_inserter.insert(SvCall(offset + 3*l - 1, offset + 2*l - 1, 0, 0, True, 1000)) # c
    sv_inserter.insert(SvCall(offset + 4*l, offset + l, 0, 0, True, 1000)) # d

    sv_inserter.close(pool)
    return get_inserter.cpp_module.id

def inversion_in_inversion_2(db_conn, dataset_name, l, offset):
    JumpRunTable(db_conn)
    SvCallerRunTable(db_conn)
    get_inserter = GetCallInserter(ParameterSetManager(), db_conn, "inversion_in_inversion_2",
                                   "the sv's that were simulated", -1)
    pool = PoolContainer(1, dataset_name)
    sv_inserter = get_inserter.execute(pool)

    sv_inserter.insert(SvCall(offset + l - 1, offset + 2*l, 0, 0, False, 1000)) # a
    sv_inserter.insert(SvCall(offset + 3*l - 1, offset + 2*l - 1, 0, 0, True, 1000)) # b
    sv_inserter.insert(SvCall(offset + 3*l, offset + l, 0, 0, True, 1000)) # c

    sv_inserter.close(pool)
    return get_inserter.cpp_module.id

def insertion_in_inversion(db_conn, dataset_name, l, offset):
    JumpRunTable(db_conn)
    SvCallerRunTable(db_conn)
    get_inserter = GetCallInserter(ParameterSetManager(), db_conn, "insertion in inversion",
                                   "the sv's that were simulated", -1)
    pool = PoolContainer(1, dataset_name)
    sv_inserter = get_inserter.execute(pool)

    sv_inserter.insert(SvCall(offset + l - 1, offset + 3*l - 1, 0, 0, True, 1000)) # a

    insertion = SvCall(offset + 2*l, offset + 2*l - 1, 0, 0, False, 1000)
    insertion.inserted_sequence = NucSeq(random_nuc_seq(l))
    sv_inserter.insert(insertion) # b

    sv_inserter.insert(SvCall(offset + 3*l, offset + l, 0, 0, True, 1000)) # c

    sv_inserter.close(pool)
    return get_inserter.cpp_module.id

def inversion(db_conn, dataset_name, l, offset):
    JumpRunTable(db_conn)
    SvCallerRunTable(db_conn)
    get_inserter = GetCallInserter(ParameterSetManager(), db_conn, "inversion",
                                   "the sv's that were simulated", -1)
    pool = PoolContainer(1, dataset_name)
    sv_inserter = get_inserter.execute(pool)

    sv_inserter.insert(SvCall(offset + l - 1, offset + 2*l - 1, 0, 0, True, 1000)) # a

    sv_inserter.insert(SvCall(offset + 2*l, offset + l, 0, 0, True, 1000)) # b

    sv_inserter.close(pool)
    return get_inserter.cpp_module.id

db_name = "perfect_alignment_caller_fail"
l = 1000
coverage = 100
read_size = l*3-1
callers = [
    (sniffles, "sniffles", sniffles_interpreter),
    #(pbSv, "pbSv", pb_sv_interpreter),
]
genome = genome_dir + "/GRCh38.p12"

if __name__ == "__main__":
    if not os.path.exists( db_prefix + db_name ):
        os.mkdir(db_prefix + db_name)
    with open(db_prefix + db_name + "/info.json", "w") as json_file:
        json_file.write("{\n")
        json_file.write("\"reference_path\":\""+genome+"\"\n")
        json_file.write("}\n")


    db_conn = DbConn({"SCHEMA": {"NAME": db_name}}) # , "FLAGS": ["DROP_ON_CLOSURE"]


    reference = Pack()
    reference.load(genome + "/ma/genome")

    chr1_len = reference.contigLengths()[0]
    #section_size = l*10
    #section_size = l*5
    #section_size = l*4
    section_size = l*3
    ref_section = "N"
    while 'n' in ref_section or 'N' in ref_section:
        offset = 241244918#random.randrange(1000, chr1_len - section_size)
        ref_section = str(reference.extract_from_to(offset, offset+section_size))

    #run_id = four_nested_svs_calls(db_conn, db_name, l, offset)
    #run_id = inversion_in_inversion(db_conn, db_name, l, offset)
    #run_id = inversion_in_inversion_2(db_conn, db_name, l, offset)
    #run_id = insertion_in_inversion(db_conn, db_name, l, offset)
    run_id = inversion(db_conn, db_name, l, offset)

    call_table = SvCallTable(db_conn)
    jump_table = SvJumpTable(db_conn) # initialize jump table
    seeds, inserts = call_table.calls_to_seeds(reference, run_id)

    if True:
        seed_printer = SeedPrinter(ParameterSetManager(), "call seed", x_range=(offset, offset+section_size),
                                   y_range=(offset, offset+section_size), do_print=False)
        seed_printer.execute(seeds)

    num_reads = (coverage * section_size) // read_size
    alignments_list = alignments_from_db(call_table, reference, run_id, read_size, num_reads,
                                         offset, offset + section_size)

    if False:
        seed_printer = SeedPrinter(ParameterSetManager(), "alignment seed", x_range=(offset, offset+section_size),
                                y_range=(0, read_size), do_print=False)
        for read, alignments in alignments_list[:1]:
            alignment_seeds = Seeds()
            for alignment in alignments:
                print(alignment.cigarString(reference, read_size))
                alignment_seeds.extend(alignment.to_seeds(reference) )
            seed_printer.execute(alignment_seeds)

    file_name = db_name + "-" + str(read_size) + "-" + str(coverage) + "-" + str(l)
    sam_file_name = sam_folder + file_name
    if True: # write alignments myself
        #alignment_to_file(alignments_list, sam_folder + "us", reference)
        alignment_to_file(alignments_list, sam_file_name, reference)
    if False: # use ngmlr
        read_to_file(alignments_list, fasta_folder + file_name)
        json_dict = { "reference_path": genome }
        read_set = {
            "fasta_file": fasta_folder + file_name + ".fasta",
            "name": "perfect_alignments_caller_fail",
            "technology": "pb"
        }
        ngmlr(read_set, sam_folder + "ngmlr.sam", json_dict)
        #ngmlr(read_set, sam_file_name + "-2.sam", json_dict)
        #mm2(read_set, sam_file_name + ".sam", json_dict)
        sam_to_bam(sam_file_name)

    with open(global_prefix + "/vcf_errors.log", "w") as error_file:
        for caller, name, interpreter in callers:
            vcf_file_path = vcf_folder + file_name + "-" + name + ".vcf"
            caller(sam_file_name + ".sorted.bam", vcf_file_path, genome)
            if not os.path.exists( vcf_file_path ):
                print("caller did not create calls: ", name)
                continue
            vcf_to_db(name, "desc", db_name, vcf_file_path, reference, interpreter, error_file)