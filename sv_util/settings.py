

main_data_folder = "/MAdata"

human_genome_dir = main_data_folder + "/genome/human/GRCh38.p12"

sv_hidden_to_aligners_data_dir = main_data_folder + "/sv_caller_analysis/svs_hidden_to_aligners"
sv_hidden_err_distrib = sv_hidden_to_aligners_data_dir + "/HG002.m84011_220902_175841_s1.GRCh38.sampled_errors.json"

ambiguities_of_atomic_sv_data_dir = main_data_folder + "/sv_caller_analysis/ambiguities_of_atomic_sv/"
sam_folder = ambiguities_of_atomic_sv_data_dir + "sam/"
fasta_folder = ambiguities_of_atomic_sv_data_dir + "fasta/"
vcf_folder = ambiguities_of_atomic_sv_data_dir + "vcf/"

home = "/usr/home/markus/workspace/"
sniffles_path = home + "Sniffles/bin/sniffles-core/sniffles"
delly_path = home + "delly/delly"
bcf_tools_path = home + "bcftools/bcftools-1.9/bcftools"
survivor_path = home + "SURVIVOR/Debug/SURVIVOR"
survivor_error_profile_path = home + "SURVIVOR/HG002_PacBio_CCS_10kb_error_profile_mm2.txt"
survivor_error_profile_path_oxf_nano = home + "SURVIVOR/NA12878_nano_error_profile_bwa.txt"
dwgsim_path = home + "DWGSIM/dwgsim"
sam_tools_pref = home + "samtools/samtools "
minimap_path = home + "minimap2/minimap2"
ngmrl_path = home + "ngmlr/bin/ngmlr-0.2.7/ngmlr"
vg_path = home + "vg"
graph_aligner_path = "/usr/home/markus/miniconda3/envs/ma/bin/GraphAligner"
graph_aligner2_path = home + "GraphAligner/bin/GraphAligner"
gridss_path = home + "gridss/gridss -j ~/workspace/gridss/gridss-2.12.2-gridss-jar-with-dependencies.jar "

reconstructed_query_genome_path = main_data_folder + "/genome/reconstructed/yeast/UFRJ50816"

tmp_ksw_file_prefix = main_data_folder + "/tmp/.CIGARMemoryManager"

accuracy_recall_data_dir = main_data_folder + "/sv_caller_analysis/yeast_analysis"
gridss_data_dir = main_data_folder + "/sv_caller_analysis/yeast_analysis/caller"

# Warning: enabling the NW alignments will cause the machine to use large
# amounts of Disk & RAM space 
run_ksw = True
ksw_file_system_min_gb_size = 8

show_plots = False
save_plots = True