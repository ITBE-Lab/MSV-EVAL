import os
from sv_util.settings import *
from subprocess import run

TIMEOUT=60*30 # 30 minutes time for each alignment

def bwa(read_set, sam_file_path, json_dict):
    run(["sv_util/bwa.sh", read_set["name"], json_dict["reference_path"] + "/bwa/genome", read_set["fasta_file"], 
        read_set["fasta_file_mate"], sam_file_path], timeout=TIMEOUT)

def bwa_single(read_set, sam_file_path, json_dict):
    run(["sv_util/bwa_single.sh", read_set["name"], json_dict["reference_path"] + "/bwa/genome", 
         read_set["fasta_file"], 
         sam_file_path], timeout=TIMEOUT)

def bowtie(read_set, sam_file_path, json_dict):
    run(["sv_util/bowtie.sh", read_set["name"], json_dict["reference_path"] + "/bowtie/genome.fna", 
         read_set["fasta_file"], read_set["fasta_file_mate"], 
         sam_file_path], timeout=TIMEOUT)

def mm2(read_set, sam_file_path, json_dict, extra="", presetting_override=None):
    presetting = None # noop
    if read_set["technology"] == "pb":
        presetting = "map-pb"
    if read_set["technology"] == "ont":
        presetting = "map-ont"
    if not presetting_override is None:
        presetting = presetting_override
    index_str = json_dict["reference_path"] + "/minimap/genome." + presetting + ".mmi"
    
    run(["sv_util/mm2.sh", minimap_path, presetting, extra, read_set["name"], index_str, 
         read_set["fasta_file"], sam_file_path], 
        timeout=TIMEOUT)

def pbmm2(read_set, sam_file_path, json_dict):
    run(["sv_util/pbmm2.sh", "pbmm2", json_dict["reference_path"] + "/minimap/genome.map-pb.mmi",
         read_set["fasta_file"], sam_file_path + ".bam"], timeout=TIMEOUT)
    s = sam_tools_pref + "view -h -o " + sam_file_path + " " + sam_file_path + ".bam"
    run(s, shell=True)

def sw(read_set, sam_file_path, json_dict):
    run(["sv_util/sw.sh", json_dict["reference_path"] + "/fasta/genome.fna", read_set["fasta_file"], sam_file_path], 
        timeout=TIMEOUT)

def ngmlr(read_set, sam_file_path, json_dict):
    presetting = None # noop
    if read_set["technology"] == "pb":
        presetting = "pacbio"
    if read_set["technology"] == "ont":
        presetting = "ont"
    index_str = json_dict["reference_path"] + "/ngmlr/genome.fna"
    run(["sv_util/ngmlr.sh", ngmrl_path, index_str, read_set["fasta_file"], 
         read_set["name"], presetting, sam_file_path], timeout=TIMEOUT)

def blasr(read_set, sam_file_path, json_dict):
    #print(s)
    run(["sv_util/blasr.sh", read_set["fasta_file"], json_dict["reference_path"], sam_file_path + ".no_qstring"], 
        timeout=TIMEOUT)
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

def graph_aligner(read_set, sam_file_path, json_dict):
    #print(s)
    run(["sv_util/graph_aligner.sh", graph_aligner_path, json_dict["reference_path"] + "/vg/genome.vg", 
         read_set["fasta_file"], sam_file_path + ".gam"], 
        timeout=TIMEOUT)            
    run(["sv_util/vg.sh", vg_path, json_dict["reference_path"] + "/vg/genome.vg", 
         sam_file_path + ".gam", sam_file_path + ".bam"], 
        timeout=TIMEOUT)
    #print(s)
    s = sam_tools_pref + "view -h -o " + sam_file_path + " " + sam_file_path + ".bam 2>/dev/null"
    #print(s)
    run(s, shell=True, timeout=TIMEOUT)

def graph_aligner_2(read_set, sam_file_path, json_dict):
    #print(s)
    run(["sv_util/graph_aligner.sh", graph_aligner2_path, json_dict["reference_path"] + "/vg/genome.vg", 
         read_set["fasta_file"], sam_file_path + ".gam"], 
        timeout=TIMEOUT)            
    run(["sv_util/vg.sh", vg_path, json_dict["reference_path"] + "/vg/genome.vg", 
         sam_file_path + ".gam", sam_file_path + ".bam"], 
        timeout=TIMEOUT)
    #print(s)
    s = sam_tools_pref + "view -h -o " + sam_file_path + " " + sam_file_path + ".bam 2>/dev/null"
    #print(s)
    run(s, shell=True, timeout=TIMEOUT)

def sam_to_bam(sam_file_path):
    dev_null = "" #" 2> /dev/null "
    # create sorted and indexed bam files
    to_bam_cmd = sam_tools_pref + "view -Sb " + sam_file_path + ".sam > " + sam_file_path + ".bam"
    os.system(to_bam_cmd)
    sort_cmd = sam_tools_pref + "sort -@ 7 -m 1G " + sam_file_path + ".bam > " \
                + sam_file_path + ".sorted.bam"
    os.system(sort_cmd + dev_null)
    index_cmd = sam_tools_pref + "index " + sam_file_path + ".sorted.bam > " \
                + sam_file_path + ".sorted.bam.bai"
    os.system(index_cmd)