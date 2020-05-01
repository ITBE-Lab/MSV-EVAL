from MSV import *
from MA import *
import random
from perfect_alignments_caller_fail.alignments_from_db import *
from bokeh.io.output import reset_output

global_prefix = "/MAdata/sv_lost_during_calling/"
sam_folder = global_prefix + "sam/"
vcf_folder = global_prefix + "vcf/"
genome_dir = "/MAdata/genome/human"

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
    get_inserter = GetCallInserter(ParameterSetManager(), db_conn, "simulated_sv",
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

db_name = "perfect_alignment_caller_fail"
if __name__ == "__main__":
    db_conn = DbConn({"SCHEMA": {"NAME": db_name, "FLAGS": ["DROP_ON_CLOSURE"]}})


    reference = Pack()
    reference.load(genome_dir + "/GRCh38.p12" + "/ma/genome")

    chr1_len = reference.contigLengths()[0]
    l = 100
    ref_section = "N"
    while 'n' in ref_section or 'N' in ref_section:
        offset = random.randrange(1000, chr1_len - l*10)
        ref_section = str(reference.extract_from_to(offset, offset+10*l))

    run_id = four_nested_svs_calls(db_conn, db_name, l, offset)

    call_table = SvCallTable(db_conn)
    seeds, inserts = call_table.calls_to_seeds(reference, run_id)

    seed_printer = SeedPrinter(ParameterSetManager(), "call seed", x_range=(offset, offset+10*l),
                               y_range=(offset, offset+10*l), do_print=False)
    seed_printer.execute(seeds)

    read_size = l*7
    alignments_list = alignments_from_db(call_table, reference, run_id, read_size, 3, offset, offset + 10*l)

    seed_printer2 = SeedPrinter(ParameterSetManager(), "alignment seed", x_range=(offset, offset+10*l),
                               y_range=(0, read_size), do_print=False)
    for read, alignments in alignments_list[:1]:
        alignment_seeds = Seeds()
        for alignment in alignments:
            print(alignment.cigarString(reference, read_size))
            alignment_seeds.extend(alignment.to_seeds(reference) )
        seed_printer2.execute(alignment_seeds)