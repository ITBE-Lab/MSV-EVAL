

main_data_folder = "/MAdata"

human_genome_dir = main_data_folder + "/genome/human/GRCh38.p12"


sv_hidden_to_aligners_data_dir = main_data_folder + "/sv_caller_analysis/svs_hidden_to_aligners"

ambiguities_of_atomic_sv_data_dir = main_data_folder + "/sv_caller_analysis/ambiguities_of_atomic_sv/"
sam_folder = ambiguities_of_atomic_sv_data_dir + "sam/"
fasta_folder = ambiguities_of_atomic_sv_data_dir + "fasta/"
vcf_folder = ambiguities_of_atomic_sv_data_dir + "vcf/"

sniffles_path = "~/workspace/Sniffles/bin/sniffles-core/sniffles"
delly_path = "~/workspace/delly/delly"
bcf_tools_path = "~/workspace/bcftools/bcftools-1.9/bcftools"
survivor_path = "~/workspace/SURVIVOR/Debug/SURVIVOR"
survivor_error_profile_path = "~/workspace/SURVIVOR/HG002_PacBio_CCS_10kb_error_profile_mm2.txt"
dwgsim_path = "~/workspace/DWGSIM/dwgsim"
sam_tools_pref = "~/workspace/samtools/samtools "
minimap_path = "~/workspace/minimap2/minimap2"
ngmrl_path = "~/workspace/ngmlr/bin/ngmlr/ngmlr"

reconstructed_query_genome_path = main_data_folder + "/genome/reconstructed/yeast/UFRJ50816"

tmp_ksw_file_prefix = main_data_folder + "/tmp/.CIGARMemoryManager"

accuracy_recall_data_dir = main_data_folder + "/sv_caller_analysis/yeast_analysis"

ksw_file_system_min_gb_size = 8