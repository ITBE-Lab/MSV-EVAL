from MSV import *
from MA import *
import random
from perfect_alignments_caller_fail.alignments_from_db import *
from bokeh.io.output import reset_output
from caller_analysis.os_sv_callers import *
from caller_analysis.vcf_interpreters import *
import os
from sv_util.os_aligners import *
from perfect_alignments_caller_fail.compute_errors import *

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

    sv_inserter.insert(SvCall(offset + l - 1, offset + 9*l - 1, 0, 0, True, False, 1000)) # a
    sv_inserter.insert(SvCall(offset + 4*l - 1, offset + 8*l, 0, 0, True, True, 1000)) # b

    insertion = SvCall(offset + 3*l - 1, offset + 3*l, 0, 0, True, True, 1000)
    insertion.inserted_sequence = NucSeq(random_nuc_seq(l))
    sv_inserter.insert(insertion) # c

    sv_inserter.insert(SvCall(offset + 2*l, offset + 5*l - 1, 0, 0, False, False, 1000)) # d
    sv_inserter.insert(SvCall(offset + 4*l, offset + 8*l - 1, 0, 0, False, False, 1000)) # e
    sv_inserter.insert(SvCall(offset + 6*l - 1, offset + 7*l, 0, 0, True, True, 1000)) # f
    sv_inserter.insert(SvCall(offset + 2*l - 1, offset + 5*l, 0, 0, True, True, 1000)) # g
    sv_inserter.insert(SvCall(offset + l, offset + 9*l, 0, 0, False, True, 1000)) # h

    sv_inserter.close(pool)
    return get_inserter.cpp_module.id

def inversion_in_inversion(db_conn, dataset_name, l, offset):
    JumpRunTable(db_conn)
    SvCallerRunTable(db_conn)
    get_inserter = GetCallInserter(ParameterSetManager(), db_conn, "inversion_in_inversion",
                                   "the sv's that were simulated", -1)
    pool = PoolContainer(1, dataset_name)
    sv_inserter = get_inserter.execute(pool)

    sv_inserter.insert(SvCall(offset + l - 1, offset + 4*l - 1, 0, 0, True, False, 1000)) # a
    sv_inserter.insert(SvCall(offset + 2*l, offset + 3*l, 0, 0, False, True, 1000)) # b
    sv_inserter.insert(SvCall(offset + 2*l - 1, offset + 3*l - 1, 0, 0, True, False, 1000)) # c
    sv_inserter.insert(SvCall(offset + l, offset + 4*l, 0, 0, False, True, 1000)) # d

    sv_inserter.close(pool)
    return get_inserter.cpp_module.id

def inversion_in_inversion_2(db_conn, dataset_name, l, offset):
    JumpRunTable(db_conn)
    SvCallerRunTable(db_conn)
    get_inserter = GetCallInserter(ParameterSetManager(), db_conn, "inversion_in_inversion_2",
                                   "the sv's that were simulated", -1)
    pool = PoolContainer(1, dataset_name)
    sv_inserter = get_inserter.execute(pool)

    sv_inserter.insert(SvCall(offset + l - 1, offset + 2*l, 0, 0, True, True, 1000)) # a
    sv_inserter.insert(SvCall(offset + 2*l - 1, offset + 3*l - 1, 0, 0, True, False, 1000)) # b
    sv_inserter.insert(SvCall(offset + l, offset + 3*l, 0, 0, False, True, 1000)) # c

    sv_inserter.close(pool)
    return get_inserter.cpp_module.id

def insertion_in_inversion(db_conn, dataset_name, l, offset):
    JumpRunTable(db_conn)
    SvCallerRunTable(db_conn)
    get_inserter = GetCallInserter(ParameterSetManager(), db_conn, "insertion in inversion",
                                   "the sv's that were simulated", -1)
    pool = PoolContainer(1, dataset_name)
    sv_inserter = get_inserter.execute(pool)

    sv_inserter.insert(SvCall(offset + l - 1, offset + 3*l - 1, 0, 0, True, False, 1000)) # a

    insertion = SvCall(offset + 2*l-1, offset + 2*l, 0, 0, True, True, 1000)
    insertion.inserted_sequence = NucSeq(random_nuc_seq(l))
    sv_inserter.insert(insertion) # b

    sv_inserter.insert(SvCall(offset + l, offset + 3*l, 0, 0, False, True, 1000)) # c

    sv_inserter.close(pool)
    return get_inserter.cpp_module.id

def inversion(db_conn, dataset_name, l, offset):
    JumpRunTable(db_conn)
    SvCallerRunTable(db_conn)
    get_inserter = GetCallInserter(ParameterSetManager(), db_conn, "inversion",
                                   "the sv's that were simulated", -1)
    pool = PoolContainer(1, dataset_name)
    sv_inserter = get_inserter.execute(pool)

    sv_inserter.insert(SvCall(offset + l - 1, offset + 2*l - 1, 0, 0, True, False, 1000)) # a

    sv_inserter.insert(SvCall(offset + l, offset + 2*l, 0, 0, False, True, 1000)) # b

    sv_inserter.close(pool)
    return get_inserter.cpp_module.id

def inversion_in_translocation(db_conn, dataset_name, l, offset):
    JumpRunTable(db_conn)
    SvCallerRunTable(db_conn)
    get_inserter = GetCallInserter(ParameterSetManager(), db_conn, "touching inversion in translocation",
                                   "the sv's that were simulated", -1)
    pool = PoolContainer(1, dataset_name)
    sv_inserter = get_inserter.execute(pool)

    sv_inserter.insert(SvCall(offset + l - 1, offset + 3*l, 0, 0, True, True, 1000)) # a
    sv_inserter.insert(SvCall(offset + 3*l-1, offset + 4*l-1, 0, 0, True, False, 1000)) # b
    sv_inserter.insert(SvCall(offset + l, offset + 2*l, 0, 0, False, True, 1000)) # c
    sv_inserter.insert(SvCall(offset + 2*l-1, offset + 4*l, 0, 0, True, True, 1000)) # d

    sv_inserter.close(pool)
    return get_inserter.cpp_module.id

def proper_inversion_in_translocation(db_conn, dataset_name, l, offset):
    JumpRunTable(db_conn)
    SvCallerRunTable(db_conn)
    get_inserter = GetCallInserter(ParameterSetManager(), db_conn, "proper inversion in translocation",
                                   "the sv's that were simulated", -1)
    pool = PoolContainer(1, dataset_name)
    sv_inserter = get_inserter.execute(pool)

    sv_inserter.insert(SvCall(offset + l - 1, offset + 5*l, 0, 0, True, True, 1000)) # a
    sv_inserter.insert(SvCall(offset + 2*l, offset + 6*l-1, 0, 0, False, False, 1000)) # b
    sv_inserter.insert(SvCall(offset + 3*l-1, offset + 4*l-1, 0, 0, True, False, 1000)) # c
    sv_inserter.insert(SvCall(offset + 3*l, offset + 4*l, 0, 0, False, True, 1000)) # d
    sv_inserter.insert(SvCall(offset + l, offset + 5*l-1, 0, 0, False, False, 1000)) # e
    sv_inserter.insert(SvCall(offset + 2*l-1, offset + 6*l, 0, 0, True, True, 1000)) # f

    sv_inserter.close(pool)
    return get_inserter.cpp_module.id

db_name = "perfect_alignment_caller_fail"
l = 1000
coverage = 100
read_size = l*5-1
callers = [
    (sniffles, "sniffles", sniffles_interpreter, "single"),
    #(pbSv, "pbSv", pb_sv_interpreter),
]
genome = genome_dir + "/GRCh38.p12"

def run_msv(pack, alignments_list):
    params = ParameterSetManager()

    jumps_from_seeds = SvJumpsFromSeeds(params, pack)
    get_jump_inserter = GetJumpInserter(params, DbConn(db_name), "MS-SV", "jumps from perfect seeds")
    jump_inserter_module = JumpInserterModule(params)
    pool = PoolContainer(1, db_name)
    jump_inserter = get_jump_inserter.execute(pool)

    for query, alignments in alignments_list:
        alignment_seeds = Seeds()
        for alignment in alignments:
            alignment_seeds.extend(alignment.to_seeds(reference))
        jumps = jumps_from_seeds.cpp_module.compute_jumps(alignment_seeds, query, pack)

        if False:
            #if alignment_seeds[0].on_forward_strand:
            #    continue
            seed_printer = SeedPrinter(ParameterSetManager(), "msv seeds", do_print=True)
            seed_printer.execute(alignment_seeds, None,
                                 jumps_from_seeds.cpp_module.execute_helper_no_reseed(alignment_seeds, pack, query))
            exit()

        jump_inserter_module.execute(jump_inserter, pool, jumps, query)
    jump_inserter.close(pool)
    jump_id = get_jump_inserter.cpp_module.id
    return sweep_sv_jumps(params, db_name, jump_id, "MA_SV", "", [-1], pack, silent=True)


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

    svs = [
        (l*10, four_nested_svs_calls),
        (l*5, inversion_in_inversion),
        (l*4, inversion_in_inversion_2),
        (l*4, insertion_in_inversion),
        (l*3, inversion),
        (l*5, inversion_in_translocation),
        (l*7, proper_inversion_in_translocation),
    ]

    sets = []
    for section_size, sc_func in svs:
        ref_section = "N"
        while 'n' in ref_section or 'N' in ref_section:
            offset = random.randrange(1000, chr1_len - section_size)
            ref_section = str(reference.extract_from_to(offset, offset+section_size))

        run_id = sc_func(db_conn, db_name, l, offset)

        call_table = SvCallTable(db_conn)
        jump_table = SvJumpTable(db_conn) # initialize jump table
        seeds, inserts = call_table.calls_to_seeds(reference, run_id)

        if True:
            seed_printer = SeedPrinter(ParameterSetManager(), "call seed", x_range=(offset, offset+section_size),
                                    y_range=(offset, offset+section_size), do_print=False)
            seed_printer.execute(seeds)
            exit()

        num_reads = (coverage * section_size) // read_size
        alignments_list = alignments_from_db(call_table, reference, run_id, read_size, num_reads,
                                            offset, offset + section_size)

        if False:
            seed_printer = SeedPrinter(ParameterSetManager(), "alignment seed", x_range=(offset, offset+section_size),
                                    y_range=(0, read_size), do_print=False)
            for read, alignments in alignments_list[:1]:
                alignment_seeds = Seeds()
                for alignment in alignments:
                    print(alignment.cigarString(reference, read_size, False))
                    alignment_seeds.extend( alignment.to_seeds(reference) )
                for seed in seeds:
                    if seed.start >= offset:
                        seed.start -= offset
                seed_printer.execute(alignment_seeds, seeds)
                exit()

        file_name = db_name + "-" + str(read_size) + "-" + str(coverage) + "-" + str(l)
        sam_file_name = sam_folder + file_name
        if False: # write alignments myself
            alignment_to_file(alignments_list[:1], sam_folder + "us", reference)
        if False: # use ngmlr
            read_to_file(alignments_list[:1], fasta_folder + file_name)
            json_dict = { "reference_path": genome }
            read_set = {
                "fasta_file": fasta_folder + file_name + ".fasta",
                "name": "perfect_alignments_caller_fail",
                "technology": "pb"
            }
            ngmlr(read_set, sam_folder + "ngmlr.sam", json_dict)
            #ngmlr(read_set, sam_file_name + ".sam", json_dict)
            #mm2(read_set, sam_file_name + ".sam", json_dict)
            #sam_to_bam(sam_file_name)
        alignment_to_file(alignments_list, sam_file_name, reference)

        caller_ids = [run_msv(reference, alignments_list)]
        read_types = ["single"]
        # other callers
        with open(global_prefix + "/vcf_errors.log", "w") as error_file:
            for caller, name, interpreter, read_type in callers:
                vcf_file_path = vcf_folder + file_name + "-" + name + ".vcf"
                caller(sam_file_name + ".sorted.bam", vcf_file_path, genome)
                if not os.path.exists( vcf_file_path ):
                    print("caller did not create calls: ", name)
                    continue
                call_id = vcf_to_db(name, "desc", db_name, vcf_file_path, reference, interpreter, error_file)
                caller_ids.append(call_id)
                read_types.append(read_type)
                caller_seeds, inserts = call_table.calls_to_seeds(reference, call_id)
                if False:
                    seed_printer = SeedPrinter(ParameterSetManager(), "caller seed",
                                            x_range=(offset, offset+section_size),
                                            y_range=(offset, offset+section_size), do_print=False)
                    seed_printer.execute(caller_seeds, seeds)
        sets.append( (caller_ids, read_types, run_id) )
    render(db_conn, sets)