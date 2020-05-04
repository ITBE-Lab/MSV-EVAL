import os

def bwa(read_set, sam_file_path, json_dict):
    index_str = json_dict["reference_path"] + "/bwa/genome"
    os.system("~/workspace/bwa/bwa mem -R \"@RG\\tID:1\\tSM:" + read_set["name"] + "\" -t 32 " + index_str + " "
                + read_set["fasta_file"] + " " + read_set["fasta_file_mate"] + " > " + sam_file_path
                + " 2> /dev/null")

def bowtie(read_set, sam_file_path, json_dict):
    index_str = json_dict["reference_path"] + "/bowtie/genome.fna"
    os.system("~/workspace/bowtie2/bowtie2-2.3.3.1/bowtie2 --rg-id 1 --rg SM:" + read_set["name"] + " -p 32 -x " +
                index_str + " -1 " + read_set["fasta_file"] + " -2 " + read_set["fasta_file_mate"] + " -S " + sam_file_path
                + " 2> /dev/null")

def mm2(read_set, sam_file_path, json_dict, extra=""):
    presetting = None # noop
    if read_set["technology"] == "pb":
        presetting = "map-pb"
    if read_set["technology"] == "ont":
        presetting = "map-ont"
    index_str = json_dict["reference_path"] + "/minimap/genome." + presetting + ".mmi"
    # -c output CIGAR in PAF; -a output SAM
    s = "~/workspace/minimap2/minimap2 --MD -c -a -t 32 -x " + presetting + " " + extra + " -R \"@RG\\tID:1\\tSM:" \
        + read_set["name"] + "\" " + index_str + " " \
        + read_set["fasta_file"] + " > " + sam_file_path + " 2> /dev/null"
    #print(s)
    os.system(s)

def pbmm2(read_set, sam_file_path, json_dict):
    presetting = None # noop
    if read_set["technology"] == "pb":
        presetting = "map-pb"
    if read_set["technology"] == "ont":
        presetting = "map-ont"
    index_str = json_dict["reference_path"] + "/minimap/genome." + presetting + ".mmi"
    # -c output CIGAR in PAF; -a output SAM; -Y softclip
    os.system("~/workspace/minimap2/minimap2 --MD -Y -c -a -t 32 -x " + presetting + " -R \"@RG\\tID:1\\tSM:"
                + read_set["name"] + "\" " + index_str + " "
                + read_set["fasta_file"] + " > " + sam_file_path + " 2> /dev/null")

    #os.system("~/miniconda2/bin/pbmm2 align " + json_dict["reference_path"] + "/fasta/genome.fna " 
    #          + read_set["fasta_file"] + ".bam " + sam_file_path + 
    #          ".sorted.bam --sort --preset CCS --sample sample1 --rg '@RG\tID:movie1' ")

def ngmlr(read_set, sam_file_path, json_dict):
    presetting = None # noop
    if read_set["technology"] == "pb":
        presetting = "pacbio"
    if read_set["technology"] == "ont":
        presetting = "ont"
    index_str = json_dict["reference_path"] + "/ngmlr/genome.fna"
    # -c output CIGAR in PAF; -a output SAM
    os.system("~/workspace/ngmlr/bin/ngmlr-0.2.8/ngmlr -r " + index_str + " -q " + read_set["fasta_file"]
                + " --rg-id 1, --rg-sm " + read_set["name"]
                + " -t 32 -x " + presetting + " > " + sam_file_path ) # + " 2> /dev/null")

def blasr(read_set, sam_file_path, json_dict):
    s = "~/workspace/legacy_blasr/blasr " + read_set["fasta_file"] + " " + json_dict["reference_path"] \
                + "/blasr/genome.fasta -printSAMQV -nproc 32 -sam -out " + sam_file_path + ".no_qstring > /dev/null"
    #print(s)
    os.system(s)
    with open(sam_file_path + ".no_qstring", "r") as in_file:
        with open(sam_file_path, "w") as out_file:
            for line in in_file:
                if line[0] == "@":
                    out_file.write(line) # \n contained in line
                else:
                    columns = line[:-1].split("\t")
                    columns[10] = "~"*len(columns[9]) # overwrite qString
                    for column in columns:
                        out_file.write(column + "\t")
                    out_file.write("\n")

def sam_to_bam(sam_file_path):
    # create sorted and indexed bam files
    sam_tools_pref = "~/workspace/samtools/samtools "
    to_bam_cmd = sam_tools_pref + "view -Sb " + sam_file_path + ".sam > " + sam_file_path + ".bam"
    os.system(to_bam_cmd)
    sort_cmd = sam_tools_pref + "sort -@ 32 -m 1G " + sam_file_path + ".bam > " \
                + sam_file_path + ".sorted.bam"
    os.system(sort_cmd + " 2> /dev/null")
    index_cmd = sam_tools_pref + "index " + sam_file_path + ".sorted.bam > " \
                + sam_file_path + ".sorted.bam.bai"
    os.system(index_cmd)