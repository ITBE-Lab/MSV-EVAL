import fnmatch
import os
from MSV import *
from sv_util.settings import *

#genome_dir = global_prefix + "genome/human/GRCh38.p12/"
genome_dir = main_data_folder + "/genome/yeasts/YPS138/"
read_data_dir = main_data_folder + "/ena/" #"giab/"

def regex_match(folder, regex):
    l = []
    for file in os.listdir(folder):
        if fnmatch.fnmatch(file, regex):
            l.append(folder + str(file))
    return l


def load_reads(individual, param, read_datasets):
    seq_ids = []
    for name, f_path_vec_1, f_path_vec_2, coverage in read_datasets:
        seq_ids.append(insert_reads_path_string_vec(param, individual, name, f_path_vec_1,
                                                    f_path_vec_2, coverage=coverage))
    return seq_ids

def compute_jumps_n_calls(individual, param, seq_ids, pack, mm_index):
    jump_id = compute_sv_jumps(param, mm_index, pack, individual, seq_ids)
    sv_caller_run_id = sweep_sv_jumps(param, individual, jump_id, "MA", "", [0], pack)
    return sv_caller_run_id

def run_ma(pack, read_datasets, individual="HG002", parameter_set="SV-PacBio"):
    param = ParameterSetManager()
    param.set_selected(parameter_set)
    #param.by_name("Do Dummy Jumps").set(False) # required for real world reads
    mm_index = MinimizerIndex(param, pack.contigSeqs(), pack.contigNames())
    seq_ids = load_reads(individual, param, read_datasets)
    sv_caller_run_id = compute_jumps_n_calls(individual, param, seq_ids, pack, mm_index)
    print("caller_id", sv_caller_run_id)
    return sv_caller_run_id


def simulate_reads():
    os.system(survivor_path +
              " simreads " + main_data_folder + "/genome/yeasts/UFRJ50816/fasta/genome.fna " +
              survivor_error_profile_path + " 100 " + read_data_dir +
              "simulated/UFRJ50816/pacbio_CCS/survivor_reads.fasta")
    os.system(dwgsim_path + " -r 0 -1 250 -2 250 /MAdata/genome/yeasts/UFRJ50816/fasta/genome.fna /MAdata/ena/simulated/UFRJ50816/Illumina-250/")

if True:
    #simulate_reads()

    pack = Pack()
    pack.load(genome_dir + "ma/genome")

    ## real world reads (read types do not match with simulated ones, so obviously results vary)

    ## trim real world illumina reads
    ## java -jar ~/workspace/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar PE -threads 32 SRR4074411.1_1.fastq SRR4074411.1_2.fastq SRR4074411.1_1.paired.trimmed.fastq SRR4074411.1_1.unpaired.trimmed.fastq SRR4074411.1_2.paired.trimmed.fastq SRR4074411.1_2.unpaired.trimmed.fastq ILLUMINACLIP:adapters.fa:2:30:10 SLIDINGWINDOW:5:20 MINLEN:36
    #run_ma(pack,
    #       [("Illumina",
    #         regex_match(read_data_dir + "PRJEB7245/UFRJ50816/illumina_hiseq_2500/", "*.trimmed.fastq"),
    #         None,
    #         50),],
    #       individual="UFRJ50816",
    #       parameter_set="SV-Illumina")
    #run_ma(pack,
    #       [("PacBio",
    #         regex_match(read_data_dir + "PRJEB7245/UFRJ50816/pacBioSMRT/", "*.fasta"),
    #         None,
    #         50)],
    #       individual="UFRJ50816", 
    #       parameter_set="SV-PacBio")

    ## simulated reads
    if False:
        run_id = run_ma(pack,
                        [("SimulatedIllumina",
                            regex_match(read_data_dir + "simulated/UFRJ50816/Illumina-250/", "*.bwa.read*.fastq.gz"),
                            None,
                            100),],
                        individual="UFRJ50816",
                        parameter_set="SV-Illumina")
    else:
        run_id = 1
    SvCallTable(DbConn("UFRJ50816")).extract_small_calls(run_id, 10, "MA < 10nt",
                                                          "Illumina based MA calls smaller than 10nt")
    if False:
        run_ma(pack,
            [("SimulatedPacBio",
                regex_match(read_data_dir + "simulated/UFRJ50816/pacbio_CCS/", "*.fasta"),
                None,
                100)],
            individual="UFRJ50816", 
            parameter_set="SV-PacBio")


