#from high_confidence.vcf_interpreter_confidence import *
import fnmatch
import os
from MSV import *

global_prefix = "/MAdata/"
#genome_dir = global_prefix + "genome/human/GRCh38.p12/"
genome_dir = global_prefix + "genome/yeasts/YPS138/"
read_data_dir = global_prefix + "ena/" #"giab/"

#def load_high_confidence_calls(pack, #file_name="HG002_GRCh38_GIAB_highconf_CG-Illfb-IllsentieonHC-Ion-10XsentieonHC-SOLIDgatkHC_CHROM1-22_v.3.3.#2_highconf_triophased.vcf", individual="HG002"):
#
#    calls_id, _ = vcf_to_db("high_confidence_calls", "n/a", individual, data_dir + file_name, pack,
#                            confidence_call_interpreter)
#
#    db_conn = DbConn({"SCHEMA": {"NAME": individual}}) # , "FLAGS": ["DROP_ON_CLOSURE"]
#    call_table = SvCallTable(db_conn)
#    call_desc = CallDescTable(db_conn)
#    JumpRunTable(db_conn)
#    SequencerTable(db_conn)
#    ReadTable(db_conn)
#    SvJumpTable(db_conn)
#
#    call_table.gen_indices(calls_id)
#    call_desc.gen_index()

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
    #("Illumina", regex_match(read_data_dir + "PRJEB7245/UFRJ50816/illumina_hiseq_2500/", "*.trimmed.fastq"), None, 50),
    #("PacBio", regex_match(read_data_dir + "PRJEB7245/UFRJ50816/pacBioSMRT/", "*.fasta"), None, 50),
    
    #("SimulatedPacBio", regex_match(read_data_dir + "simulated/UFRJ50816/pacbio_CCS/", "*.fasta"), None, 100),
    ("SimulatedIllumina", regex_match(read_data_dir + "simulated/UFRJ50816/Illumina-250/", "*.bwa.read*.fastq.gz"), None, 100),
]

def uFRJ50816_filter(pack):
    chr_pos = pack.start_of_sequence("chrVII")
    return (chr_pos + 659132, chr_pos + 862941)

def load_reads(individual, param):
    seq_ids = []
    for name, f_path_vec_1, f_path_vec_2, coverage in read_datasets:
        seq_ids.append(insert_reads_path_string_vec(param, individual, name, f_path_vec_1,
                                                    f_path_vec_2, coverage=coverage))
    return seq_ids

def compute_jumps_n_calls(individual, param, seq_ids, pack, mm_index):
    jump_id = compute_sv_jumps(param, mm_index, pack, individual, seq_ids)
    sv_caller_run_id = sweep_sv_jumps(param, individual, jump_id, "MA", "", [0], pack)
    return sv_caller_run_id

def run_ma(pack, individual="HG002"):
    param = ParameterSetManager()
    param.set_selected("SV-PacBio")
    #param.set_selected("SV-Illumina")
    mm_index = MinimizerIndex(param, pack.contigSeqs(), pack.contigNames())
    seq_ids = [1]#load_reads(individual, param)
    #return
    sv_caller_run_id = compute_jumps_n_calls(individual, param, seq_ids, pack, mm_index)
    print("caller_id", sv_caller_run_id)


if __name__ == "__main__":
    pack = Pack()
    pack.load(genome_dir + "ma/genome")

    run_ma(pack, individual="UFRJ50816")

    #load_high_confidence_calls(pack, individual="UFRJ50816")


## Simulated reads:
# ./SURVIVOR simreads /MAdata/genome/yeasts/UFRJ50816/fasta/genome.fna ../HG002_PacBio_CCS_10kb_error_profile_mm2.txt 100 /MAdata/ena/simulated/UFRJ50816/pacbio_CCS/survivor_reads.fasta
# ./dwgsim -r 0 -1 250 -2 250 /MAdata/genome/yeasts/UFRJ50816/fasta/genome.fna /MAdata/ena/simulated/UFRJ50816/Illumina-250/

# java -jar ~/workspace/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 32 SRR4074411.1_1.fastq SRR4074411.1_2.fastq SRR4074411.1_1.paired.trimmed.fastq SRR4074411.1_1.unpaired.trimmed.fastq SRR4074411.1_2.paired.trimmed.fastq SRR4074411.1_2.unpaired.trimmed.fastq ILLUMINACLIP:adapters.fa:2:30:10 SLIDINGWINDOW:5:20 MINLEN:36
