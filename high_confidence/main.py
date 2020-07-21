from high_confidence.vcf_interpreter_confidence import *
import fnmatch
import os
from MSV import *

global_prefix = "/MAdata/"
genome_dir = global_prefix + "genome/human/GRCh38.p12/"
data_dir = global_prefix + "sv_caller_analysis/high_confidence_calls/"
#
# Zeus
#
if True:
    #read_data_dir = "/mnt/hdd0/giab/HG002/"
    read_data_dir = "/mnt/hdd0/ena/PRJEB19900/SAMEA2757769/"

#
# Hera
#
if False:
    read_data_dir = "/MAdata/giab/HG002/"

def load_high_confidence_calls(pack, file_name="HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf_triophased.vcf", individual="HG002"):

    calls_id, _ = vcf_to_db("high_confidence_calls", "n/a", individual, data_dir + file_name, pack,
                            confidence_call_interpreter)

    db_conn = DbConn({"SCHEMA": {"NAME": individual}}) # , "FLAGS": ["DROP_ON_CLOSURE"]
    call_table = SvCallTable(db_conn)
    call_desc = CallDescTable(db_conn)
    JumpRunTable(db_conn)
    SequencerTable(db_conn)
    ReadTable(db_conn)
    SvJumpTable(db_conn)

    call_table.gen_indices(calls_id)
    call_desc.gen_index()

def regex_match(folder, regex):
    l = []
    for file in os.listdir(folder):
        if fnmatch.fnmatch(file, regex):
            l.append(folder + str(file))
    return l


#read_datasets = [
#    ("Illumina_2x250bps", regex_match(read_data_dir+"Illumina_2x250bps/", "*_R1_*.fastq.gz"),
#     regex_match(read_data_dir+"Illumina_2x250bps/", "*_R2_*.fastq.gz")),
#    ("PacBio_CCS_10kb", regex_match(read_data_dir+"PacBio_CCS_10kb/", "m*.Q20.fastq"), None),
#]
#read_datasets = [
#    ("Illumina_2x250bps_PacBio_CCS_10kb", ["/mnt/hdd0/read_selection/HG002_illumina_pacbCss_unique_mem.fasta"], None),
#]

read_datasets = [
    ("Illumina_2x250bps", regex_match(read_data_dir, "*.fastq.gz"), None),
]

def load_reads(individual, param):
    seq_ids = []
    for name, f_path_vec_1, f_path_vec_2 in read_datasets:
        seq_ids.append(insert_reads_path_string_vec(param, individual, name, f_path_vec_1,
                                                    f_path_vec_2))
    return seq_ids

def compute_jumps_n_calls(individual, param, seq_ids, pack, fm_index):
    jump_id = 1 #compute_sv_jumps(param, fm_index, pack, individual, seq_ids)
    sv_caller_run_id = sweep_sv_jumps(param, individual, jump_id, "MA", "", [0], pack)
    return sv_caller_run_id

def run_ma(pack, fm_index, individual="HG002"):
    param = ParameterSetManager()
    param.by_name("Min Size Edge").set(min_sv_size)
    param.by_name("Maximal Ambiguity SV").set(1)
    seq_ids = [1] #load_reads(individual, param)
    sv_caller_run_id = compute_jumps_n_calls(individual, param, seq_ids, pack, fm_index)
    print("caller_id", sv_caller_run_id)


if __name__ == "__main__":
    pack = Pack()
    pack.load(genome_dir + "ma/genome")

    fm_index = FMIndex()
    fm_index.load(genome_dir + "ma/genome")

    run_ma(pack, fm_index, individual="SAMEA2757769")

    load_high_confidence_calls(pack, individual="SAMEA2757769")
