from high_confidence.vcf_interpreter_confidence import *
import fnmatch
import os

global_prefix = "/MAdata/"
genome_dir = global_prefix + "genome/human/GRCh38.p12/"
db_prefix = global_prefix + "sv_datasets/"
data_dir = global_prefix + "sv_caller_analysis/high_confidence_calls/"
read_data_dir = "/mnt/hdd0/giab/HG002/"

def load_high_confidence_calls(file_name="HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.2_highconf_triophased.vcf", individual="HG002"):
    pack = Pack()
    pack.load(genome_dir + "ma/genome")

    calls_id, _ = vcf_to_db("high_confidence_calls", "n/a", individual, data_dir + file_name, pack,
                            confidence_call_interpreter)

    db_conn = DbConn({"SCHEMA": {"NAME": individual}}) # , "FLAGS": ["DROP_ON_CLOSURE"]
    JumpRunTable(db_conn)
    call_table = SvCallTable(db_conn)
    call_desc = CallDescTable(db_conn)
    SvJumpTable(db_conn)

    call_table.gen_indices(calls_id)
    call_desc.gen_index()

def regex_match(folder, regex):
    l = []
    for file in os.listdir(folder):
        if fnmatch.fnmatch(file, regex):
            l.append(folder + str(file))
    return l


read_datasets = [
    ("Illumina_2x250bps", regex_match(read_data_dir+"Illumina_2x250bps/", "*_R1_*.fastq.gz"),
     regex_match(read_data_dir+"Illumina_2x250bps/", "*_R2_*.fastq.gz")),
    ("PacBio_CCS_10kb", regex_match(read_data_dir+"PacBio_CCS_10kb/", "m*.Q20.fastq"), None),
]

def compute_calls_ma(individual="HG002"):
    param = ParameterSetManager()
    param.by_name("Min Size Edge").set(min_sv_size)
    seq_ids = []
    for name, f_path_vec_1, f_path_vec_2 in read_datasets:
        seq_ids.append(insert_reads_path_string_vec(param, individual, name, f_path_vec_1,
                                                    f_path_vec_2))

    for seq_id in seq_ids:
        print(seq_ids)

if __name__ == "__main__":
    compute_calls_ma()
    #load_high_confidence_calls()