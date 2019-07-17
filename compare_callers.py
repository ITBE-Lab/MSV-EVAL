from MA import *
import math
import os
import json
import compute_sv_jumps
import sweep_sv_jumps
import traceback
from pathlib import Path

def create_alignments_if_necessary(dataset_name, json_dict, db, pack, fm_index, recompute_jumps=False):
    def bwa(read_set, sam_file_path):
        index_str = json_dict["reference_path"] + "/bwa/genome"
        os.system("~/workspace/bwa/bwa mem -R \"@RG\\tID:1\\tSM:" + read_set["name"] + "\" -t 32 " + index_str + " "
                  + read_set["fasta_file"] + " " + read_set["fasta_file_mate"] + " > " + sam_file_path
                  + " 2> /dev/null")

    def bowtie(read_set, sam_file_path):
        index_str = json_dict["reference_path"] + "/bowtie/genome.fna"
        os.system("~/workspace/bowtie2/bowtie2-2.3.3.1/bowtie2 --rg-id 1 --rg SM:" + read_set["name"] + " -p 32 -x " +
                  index_str + " -1 " + read_set["fasta_file"] + " -2 " + read_set["fasta_file_mate"] + " -S " + sam_file_path
                  + " 2> /dev/null")

    def mm2(read_set, sam_file_path):
        presetting = None # noop
        if read_set["technology"] == "pb":
            presetting = "map-pb"
        if read_set["technology"] == "ont":
            presetting = "map-ont"
        index_str = json_dict["reference_path"] + "/minimap/genome." + presetting + ".mmi"
        # -c output CIGAR in PAF; -a output SAM
        os.system("~/workspace/minimap2/minimap2 --MD -c -a -t 32 -x " + presetting + " -R \"@RG\\tID:1\\tSM:"
                  + read_set["name"] + "\" " + index_str + " "
                  + read_set["fasta_file"] + " > " + sam_file_path + " 2> /dev/null")

    def pbmm2(read_set, sam_file_path):
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

    def ngmlr(read_set, sam_file_path):
        presetting = None # noop
        if read_set["technology"] == "pb":
            presetting = "pacbio"
        if read_set["technology"] == "ont":
            presetting = "ont"
        index_str = json_dict["reference_path"] + "/ngmlr/genome.fna"
        # -c output CIGAR in PAF; -a output SAM
        os.system("~/workspace/ngmlr/bin/ngmlr-0.2.8/ngmlr -r " + index_str + " -q " + read_set["fasta_file"]
                    + " --rg-id 1, --rg-sm " + read_set["name"]
                    + " -t 32 -x " + presetting + " > " + sam_file_path + " 2> /dev/null")

    def blasr(read_set, sam_file_path):
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


    alignment_calls = {
        "create_reads_survivor": [mm2, ngmlr, pbmm2], #blasr @note is not required at the moment
        "create_illumina_reads_dwgsim": [bwa, bowtie]
    }

    for dataset in json_dict["datasets"]:
        for read_set in dataset["create_reads_funcs"]:
            if not "alignments" in read_set:
                read_set["alignments"] = []
            if not "jump_id" in read_set or recompute_jumps:
                print("computing jumps for MA-SV on", read_set["name"])
                params = ParameterSetManager()
                if read_set["func_name"] == "create_illumina_reads_dwgsim":
                    params.set_selected("SV-Illumina")
                elif read_set["func_name"] == "create_reads_survivor" and read_set["technology"] == "pb":
                    params.set_selected("SV-PacBio")
                elif read_set["func_name"] == "create_reads_survivor" and read_set["technology"] == "ont":
                    params.set_selected("SV-ONT")
                else:
                    print("WARNING: unknown read simulator - using default parameters for sv jumps")
                #params.by_name("Number of Threads").set(1)
                #params.by_name("Use all Processor Cores").set(False)
                read_set["jump_id"] = compute_sv_jumps.compute_sv_jumps(params, fm_index, pack, db, 
                                                                        read_set["seq_id"])
            for alignment_call in alignment_calls[read_set["func_name"]]:
                sam_file_path = "/MAdata/sv_datasets/" + dataset_name + "/alignments/" \
                            + read_set["name"] + "-" + alignment_call.__name__
                if not alignment_call.__name__ in read_set["alignments"]:
                    read_set["alignments"].append(alignment_call.__name__)
                if os.path.exists( sam_file_path + ".sam" ):
                    continue
                print("creating alignment files for", read_set["name"], alignment_call.__name__)
                alignment_call(read_set, sam_file_path + ".sam")

                #if os.path.exists( sam_file_path + ".sam.sorted.bam "):
                #    os.system("mv " + sam_file_path + ".sam.sorted.bam " + sam_file_path + ".sorted.bam ")
                #    os.system("touch " + sam_file_path + ".sam")
                if os.path.exists( sam_file_path + ".sam" ):
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
                else:
                    print("aligner did not create alignments:", alignment_call)

def vcf_parser(file_name):
    class VCFFile:
        def __init__(self, d):
            self.data = d

        def __getitem__(self, name):
            if not name in self.data:
                return []
            return self.data[name]

        def __contains__(self, name):
            return name in self.data

        def __str__(self):
            s = " {"
            for key, val in self.data.items():
                s += str(key) + ": " + str(val) + "; "
            return s + " }"

    with open(file_name, "r") as vcf_file:
        names = []
        for line in vcf_file:
            if line[-1] == "\n":
                line = line[:-1]
            if line[:2] == "##":
                pass
                #info.append(line)
            elif line[0] == "#":
                names = line[1:].split("\t")
            else:
                d = {}
                for name, field in zip(names, line.split("\t")):
                    if name == "INFO":
                        d2 = {}
                        for key_value in field.split(";"):
                            if "=" in key_value:
                                key, value = key_value.split("=")
                                d2[key] = value
                            else:
                                d2[key_value] = True
                        d[name] = VCFFile(d2)
                    else:
                        d[name] = field
                yield VCFFile(d)

def vcf_to_db(name, desc, sv_db, file_name, pack, error_file=None, create_sv_func=None):
    call_inserter = SvCallInserter(sv_db, name, desc, -1) # -1 since there are no related sv jumps...
    def find_confidence(call):
        if "coverage" in call["INFO"]:
            return int(float(call["INFO"]["coverage"]) * 10)
        if "RE" in call["INFO"]: # sniffles
            return int(float(call["INFO"]["RE"]))
        if "PE" in call["INFO"] and "SR" in call["INFO"]: # pbHoney
            return int(call["INFO"]["PE"]) + int(call["INFO"]["SR"])
        if "PE" in call["INFO"]: # pbHoney
            return int(call["INFO"]["PE"])
        if call["QUAL"] != ".": # vcf...
            return int(float(call["QUAL"]))
        return 0
    def find_std_from_std_to(call):
        std_from = None
        std_to = None
        if "STD_quant_start" in call["INFO"]:
            std_from = math.ceil(float(call["INFO"]["STD_quant_start"]))
        if "CIPOS" in call["INFO"]: # delly
            x = call["INFO"]["CIPOS"].split(",")
            std_from = math.ceil(float(x[1]) - float(x[0]))
        if "STD_quant_stop" in call["INFO"]:
            std_to = math.ceil(float(call["INFO"]["STD_quant_stop"]))
        if "CIEND" in call["INFO"]: # delly
            x = call["INFO"]["CIEND"].split(",")
            std_to = math.ceil(float(x[1]) - float(x[0]))
        return std_from, std_to
    num_calls = 0
    bnd_mate_dict = {}
    for call in vcf_parser(file_name):
        num_calls += 1
        try:
            if call["FILTER"] != "PASS": # @todo verify that this is correct
                continue
            #print("quality:", find_confidence(call))
            if call["TYPE"] == "DEL":
                #print(call)
                from_pos = int(call["START"]) + pack.start_of_sequence(call["CHROM"])
                to_pos = int(call["END"]) + pack.start_of_sequence(call["CHROM"])
                call_inserter.insert_call(SvCall(from_pos, to_pos, 1, 1, False, find_confidence(call), 1))
            elif call["ALT"] == "<DEL>" and "PRECISE" in call["INFO"]:
                #print(call)
                from_pos = int(call["POS"]) + pack.start_of_sequence(call["CHROM"])
                to_pos = int(call["INFO"]["END"]) + pack.start_of_sequence(call["INFO"]["CHR2"])
                call_inserter.insert_call(SvCall(from_pos, to_pos, 1, 1, False, find_confidence(call), 1))
            elif call["ALT"] == "<DEL>" and "IMPRECISE" in call["INFO"]:
                #print(call)
                std_from, std_to = find_std_from_std_to(call)
                from_pos = int(call["POS"]) + pack.start_of_sequence(call["CHROM"]) - int(std_from/2)
                to_pos = int(call["INFO"]["END"]) + pack.start_of_sequence(call["INFO"]["CHR2"]) - int(std_to/2)
                call_inserter.insert_call(SvCall(from_pos, to_pos, std_from, std_to, False, find_confidence(call), 1))
            elif call["ALT"] == "<DUP>" and "PRECISE" in call["INFO"]:
                #print(call)
                from_pos = int(call["POS"]) + pack.start_of_sequence(call["CHROM"])
                to_pos = int(call["INFO"]["END"]) + pack.start_of_sequence(call["INFO"]["CHR2"])
                call_inserter.insert_call(SvCall(to_pos, from_pos, 1, 1, False, find_confidence(call), 1))
            elif call["ALT"] == "<DUP>" and "IMPRECISE" in call["INFO"]:
                #print(call)
                std_from, std_to = find_std_from_std_to(call)
                from_pos = int(call["POS"]) + pack.start_of_sequence(call["CHROM"]) - int(std_from/2)
                to_pos = int(call["INFO"]["END"]) + pack.start_of_sequence(call["INFO"]["CHR2"]) - int(std_to/2)
                call_inserter.insert_call(SvCall(to_pos, from_pos, std_from, std_to, False, find_confidence(call), 1))
            elif call["ALT"] == "<DUP>": # pbsv
                #print(call)
                from_pos = int(call["POS"]) + pack.start_of_sequence(call["CHROM"])
                to_pos = from_pos + int(call["INFO"]["SVLEN"])
                call_inserter.insert_call(SvCall(to_pos, from_pos, 1, 1, False, find_confidence(call), 1))
            elif call["ALT"] == "<INS>" and "PRECISE" in call["INFO"]:
                #print(call)
                from_pos = int(call["POS"]) + pack.start_of_sequence(call["CHROM"])
                call_inserter.insert_call(SvCall(from_pos, from_pos, 1, 1, False, find_confidence(call), 1))
            elif call["TYPE"] == "INS":
                #print(call)
                from_pos = int(call["START"]) + pack.start_of_sequence(call["CHROM"])
                call_inserter.insert_call(SvCall(from_pos, from_pos, 1, 1, False, find_confidence(call), 1))
            elif call["ALT"] == "<INS>" and "IMPRECISE" in call["INFO"]:
                #print(call)
                from_start = int(call["POS"]) + pack.start_of_sequence(call["CHROM"])
                from_end = int(call["INFO"]["END"]) + pack.start_of_sequence(call["INFO"]["CHR2"])
                call_inserter.insert_call(SvCall(from_start, from_start, from_end - from_start,
                                                    from_end - from_start, False, find_confidence(call), 1))
            elif call["INFO"]["SVTYPE"] == "INS": #pbsv
                #print(call)
                from_start = int(call["POS"]) + pack.start_of_sequence(call["CHROM"])
                call_inserter.insert_call(SvCall(from_start, from_start, 1, 1, False, find_confidence(call), 1))
            elif call["ALT"] == "<INV>" and "PRECISE" in call["INFO"]:
                #print(call)
                from_pos = int(call["POS"]) + pack.start_of_sequence(call["CHROM"])
                to_pos = int(call["INFO"]["END"]) + pack.start_of_sequence(call["INFO"]["CHR2"])
                call_inserter.insert_call(SvCall(from_pos, to_pos, 1, 1, True, find_confidence(call), 1))
                call_inserter.insert_call(SvCall(to_pos, from_pos, 1, 1, True, find_confidence(call), 1))
            elif call["ALT"] == "<INV>" and "IMPRECISE" in call["INFO"]:
                #print(call)
                std_from, std_to = find_std_from_std_to(call)
                from_pos = int(call["POS"]) + pack.start_of_sequence(call["CHROM"]) - int(std_from/2)
                to_pos = int(call["INFO"]["END"]) + pack.start_of_sequence(call["INFO"]["CHR2"]) - int(std_to/2)
                call_inserter.insert_call(SvCall(from_pos, to_pos, std_from, to_pos, True, find_confidence(call), 1))
                call_inserter.insert_call(SvCall(to_pos, from_pos, from_pos, to_pos, True, find_confidence(call), 1))
            elif call["ALT"] == "<INV>": # pbsv
                #print(call)
                from_pos = int(call["POS"]) + pack.start_of_sequence(call["CHROM"])
                to_pos = from_pos + int(call["INFO"]["SVLEN"])
                call_inserter.insert_call(SvCall(from_pos, to_pos, 1, 1, True, find_confidence(call), 1))
                call_inserter.insert_call(SvCall(to_pos, from_pos, 1, 1, True, find_confidence(call), 1))
            elif call["INFO"]["SVTYPE"] == "DEL": # Manta
                from_pos = int(call["POS"]) + pack.start_of_sequence(call["CHROM"])
                to_pos = from_pos - int(call["INFO"]["SVLEN"])
                call_inserter.insert_call(SvCall(from_pos, to_pos, 1, 1, False, find_confidence(call), 1))
            elif call["INFO"]["SVTYPE"] == "DUP" and "IMPRECISE" in call["INFO"]: # Manta
                std_from, std_to = find_std_from_std_to(call)
                from_pos = int(call["POS"]) + pack.start_of_sequence(call["CHROM"]) - int(std_from/2)
                to_pos = from_pos + int(call["INFO"]["SVLEN"]) - int(std_to/2)
                call_inserter.insert_call(SvCall(to_pos, from_pos, std_from, std_to, False, find_confidence(call), 1))
            elif call["INFO"]["SVTYPE"] == "DUP": # Manta
                from_pos = int(call["POS"]) + pack.start_of_sequence(call["CHROM"])
                to_pos = from_pos + int(call["INFO"]["SVLEN"])
                call_inserter.insert_call(SvCall(to_pos, from_pos, 1, 1, False, find_confidence(call), 1))
            elif call["INFO"]["SVTYPE"] == "BND" and "IMPRECISE" in call["INFO"]: # Manta
                if call["INFO"]["MATEID"] in bnd_mate_dict:
                    mate = bnd_mate_dict[call["INFO"]["MATEID"]]
                    std_from, _ = find_std_from_std_to(mate)
                    std_to, _ = find_std_from_std_to(call)
                    from_pos = int(mate["POS"]) + pack.start_of_sequence(mate["CHROM"]) - int(std_from/2)
                    to_pos = int(call["POS"]) + pack.start_of_sequence(call["CHROM"]) - int(std_to/2)
                    save_as_inversion = create_sv_func == "sv_inversion"
                    call_inserter.insert_call(SvCall(to_pos, from_pos, std_from, std_to, save_as_inversion,
                                                     find_confidence(call), 1))
                    del bnd_mate_dict[call["INFO"]["MATEID"]]
                else:
                    bnd_mate_dict[call["ID"]] = call
            else:
                print("unrecognized sv:", call)
                error_file.write("unrecognized sv: \n")
                error_file.write(str(call))
                error_file.write("\n\n\n")
                #exit(0)
        except Exception as e:
            print("error while handeling sv:", call)
            error_file.write("error while handeling sv: \n")
            error_file.write(str(call))
            error_file.write("\n")
            error_file.write(str(e))
            error_file.write(traceback.format_exc())
            error_file.write("\n\n\n")
    for key, value in bnd_mate_dict.items():
        print("WARNING did not find mate for call", value)
        error_file.write("WARNING did not find mate for call: \n")
        error_file.write(str(value))
        error_file.write("\n\n\n")
    print("number of calls:", num_calls)
    del call_inserter # trigger deconstructor

def run_callers_if_necessary(dataset_name, json_dict, db, pack, fm_index):
    def sniffles(bam_file, vcf_file):
        os.system("~/workspace/Sniffles/bin/sniffles-core-1.0.8/sniffles -t 32 -m " + bam_file + " -v " + vcf_file
                  + " >/dev/null 2>&1")

    def pbHoney(bam_file, vcf_file):
        os.system("~/workspace/pbhoney/PBSuite_15.8.24/bin/Honey.py pie -n 32 -o " + vcf_file + ".pie.bam " + bam_file +
                  " " + json_dict["reference_path"] + "/blasr/genome.fasta")

        os.system("~/workspace/samtools/samtools sort -@ 32 -m 1G " + vcf_file + ".pie.bam > " + 
                   vcf_file + ".pie.sorted.bam")
        os.system("~/workspace/samtools/samtools index " + vcf_file + ".pie.sorted.bam > " + 
                   vcf_file + ".pie.sorted.bam.bai")

        with open(vcf_file + ".hon.tails", "w") as out_file:
            out_file.write("##fileformat=VCFv4.0\n")
            out_file.write("#CHROM\tOUTERSTART\tSTART\tINNERSTART\tINNEREND\tEND\tOUTEREND\tTYPE\tSIZE\tINFO\n")
        os.system("~/workspace/pbhoney/PBSuite_15.8.24/bin/Honey.py tails -o " + vcf_file + ".hon.tails " + 
                   vcf_file + ".pie.sorted.bam")

        with open(vcf_file + ".hon.spots.spots", "w") as out_file:
            out_file.write("##fileformat=VCFv4.0\n")
            out_file.write("#CHROM\tOUTERSTART\tSTART\tINNERSTART\tINNEREND\tEND\tOUTEREND\tTYPE\tSIZE\tINFO\n")
        #os.system("~/workspace/pbhoney/PBSuite_15.8.24/bin/Honey.py spots -n 32 -o " + vcf_file + 
        #          ".hon.spots --reference " + json_dict["reference_path"] + "/blasr/genome.fasta " + bam_file )
        os.system("~/workspace/pbhoney/PBSuite_15.8.24/bin/Honey.py spots --consensus None -n 32 -o " + vcf_file + 
                  ".hon.spots " + bam_file )

        with open(vcf_file, "w") as out_file:
            with open(vcf_file + ".hon.spots.spots", "r") as in_file:
                for line in in_file:
                    out_file.write(line)
            with open(vcf_file + ".hon.tails", "r") as in_file:
                for line in in_file:
                    if line[0] != "#":
                        out_file.write(line)


    def delly(bam_file, vcf_file):
        os.system("~/workspace/delly/delly call -g " + json_dict["reference_path"] + "/fasta/genome.fna " + bam_file
                  + " -o " + vcf_file + ".bcf ") #>/dev/null 2>&1

        os.system("~/workspace/bcftools/bcftools-1.9/bcftools view " + vcf_file + ".bcf > " + vcf_file)

    def manta(bam_file, vcf_file):
        # prepare manta
        os.system("rm -r " + vcf_file + ".manta") # clean up folder
        os.system("python2 ~/workspace/manta/manta-1.5.0.centos6_x86_64/bin/configManta.py --referenceFasta " + json_dict["reference_path"] + "/fasta/genome.fna --bam " + bam_file + " --runDir " + vcf_file + ".manta" )
        # actually run manta
        os.system("python2 " + vcf_file + ".manta/runWorkflow.py -j 32 -m local > /dev/null" )

        os.system("cp " + vcf_file + ".manta/results/variants/diploidSV.vcf.gz" + " " + vcf_file + ".gz")
        os.system("gunzip -f " + vcf_file + ".gz")

    def smoove(bam_file, vcf_file):
        docker = True
        if docker:
            bam_folder = bam_file[:bam_file.rfind("/")]
            bam_filename = bam_file[bam_file.rfind("/")+1:]
            vcf_folder = vcf_file[:vcf_file.rfind("/")]
            os.system( "rm -f " + vcf_folder + "/*.disc.*" )
            vcf_filename = vcf_file[vcf_file.rfind("/")+1:]
            s = "docker run -v " + bam_folder + ":/bam_folder/ -v " + json_dict["reference_path"] + \
                "/fasta:/genome_folder/ -v " + vcf_folder + \
                ":/vcf_folder/ -it brentp/smoove smoove call -o /vcf_folder/ --name " + vcf_filename + \
                " --noextrafilters --fasta /genome_folder/genome.fna -p 32 --genotype /bam_folder/" + bam_filename \
                + " > /dev/null"
            #print(s)
            os.system( s )
        else:
            os.system( "smoove call -o " + vcf_file + " --noextrafilters --fasta " + json_dict["reference_path"] 
                       + "/fasta/genome.fna -p 32 --genotype " + bam_file )
        os.system("gunzip -f " + vcf_file + "-smoove.genotyped.vcf.gz")

    def pbSv(bam_file, vcf_file):
        os.system("~/miniconda2/bin/pbsv discover " + bam_file + " " + vcf_file + ".svsig.gz")
        os.system("~/miniconda2/bin/pbsv call " + json_dict["reference_path"] + "/fasta/genome.fa " + 
                  vcf_file + ".svsig.gz " + vcf_file)


    sv_calls = {
        "bwa":    [delly, smoove, manta],
        "bowtie": [delly, smoove, manta],
        "mm2":    [sniffles],
        "pbmm2":  [pbSv],
        "ngmlr":  [sniffles],
        "blasr":  [] #pbHoney @note the vcf file reading seems to be broken anyways
    }

    with open("/MAdata/sv_datasets/" + dataset_name + "/vcf_errors.log", "a") as error_file:
        for dataset in json_dict["datasets"]:
            for read_set in dataset["create_reads_funcs"]:
                if not "calls" in read_set:
                    read_set["calls"] = []

                # MA-SV
                if not "MA_SV" in read_set["calls"]:
                    read_set["calls"].append("MA_SV")
                print("creating calls for", read_set["name"], "MA_SV")
                params = ParameterSetManager()
                if read_set["func_name"] == "create_illumina_reads_dwgsim":
                    params.set_selected("SV-Illumina")
                elif read_set["func_name"] == "create_reads_survivor" and read_set["technology"] == "pb":
                    params.set_selected("SV-PacBio")
                elif read_set["func_name"] == "create_reads_survivor" and read_set["technology"] == "ont":
                    params.set_selected("SV-ONT")
                else:
                    print("WARNING: unknown read simulator - using default parameters for sv jumps")
                sweep_sv_jumps.sweep_sv_jumps(params, db, read_set["jump_id"], 
                                            pack.unpacked_size_single_strand, read_set["name"] + "--" + "MA_SV",
                                            "ground_truth=" + str(dataset["ground_truth"]), [read_set["seq_id"]], 
                                            pack, fm_index)
                # other callers
                for alignment in read_set["alignments"]:
                    for sv_call in sv_calls[alignment]:
                        vcf_file_path = "/MAdata/sv_datasets/" + dataset_name + "/calls/" \
                                + read_set["name"] + "-" + alignment + "-" + sv_call.__name__ + ".vcf"
                        bam_file_path = "/MAdata/sv_datasets/" + dataset_name + "/alignments/" \
                                + read_set["name"] + "-" + alignment + ".sorted.bam"
                        if os.path.exists( vcf_file_path ):
                            print("not creating calls for", read_set["name"], alignment, sv_call.__name__)
                            continue
                        print("creating calls for", read_set["name"], alignment, sv_call.__name__)
                        sv_call(bam_file_path, vcf_file_path)
                        if not os.path.exists( vcf_file_path ):
                            print("caller did not create calls: ", read_set["name"], alignment, sv_call.__name__)
                            continue
                        vcf_to_db(read_set["name"] + "-" + alignment + "-" + sv_call.__name__,
                                 "ground_truth=" + str(dataset["ground_truth"]), db, vcf_file_path, pack, error_file,
                                 dataset["func_name"])

                    for sv_call in sv_calls[alignment]:
                        if not sv_call.__name__ in read_set["calls"]:
                            read_set["calls"].append(sv_call.__name__)


def print_columns(data):
    col_width = [max([len(data[j][i]) for j in range(len(data))]) for i in range(len(data[0]))]
    first = True
    last_row = col_width
    for row in data:
        #if not first and cat != row[0]:
        #    print("| " + "".join(" "*l + " | " for l in col_width))
        #    cat = row[0]
        print("| " + "".join(
                (word.ljust(col_width[i]) if last_row[:i+1] != row[:i+1] \
                                          else " "*col_width[i]) + " | " for i, word in enumerate(row)
              ) )
        if first:
            print("-" * (sum(col_width) + len(col_width)*3 + 1))
            first = False
        last_row = row

blur_amount = 10

def analyze_by_score(sv_db, id_a, id_b):
    num_calls_a = sv_db.get_num_calls(id_a, 0) # num calls made
    if num_calls_a == 0:
        return [], [], [], 0, 0
    min_score = sv_db.get_min_score(id_a)
    max_score = sv_db.get_max_score(id_a)
    p = min_score
    inc = (max_score - min_score) / 50
    if max_score <= min_score:
        inc = 1
    #xs = []
    xs_2 = []
    #ys = []
    ys_2 = []
    ps = []
    num_calls_b = sv_db.get_num_calls(id_b, 0) # num actual calls
    if num_calls_b == 0 or min_score == float('inf') or max_score == float('inf'):
        #print(num_calls_a, min_score, max_score)
        return [], [], [], 0, 0
    while p <= max_score:
        #print(min_score, p, max_score)
        ps.append(p)
        # how many of the sv's are detected?
        #num_overlaps_b_to_a = sv_db.get_num_overlaps_between_calls(id_b, id_a, p, 0)
        num_almost_overlaps_b_to_a = sv_db.get_num_overlaps_between_calls(id_b, id_a, p, blur_amount)

        #xs.append(num_overlaps_b_to_a/num_calls_b)
        xs_2.append(num_almost_overlaps_b_to_a/num_calls_b)

        num_calls_a = sv_db.get_num_calls(id_a, p) # num calls made

        if num_calls_a == 0:
            #ys.append(0)
            ys_2.append(0)
        else:
            #ys.append(num_overlaps_b_to_a/num_calls_a)
            ys_2.append(num_almost_overlaps_b_to_a/num_calls_a)

        p += inc

    #num_invalid_calls = sv_db.get_num_invalid_calls(id_a, 0, 0)
    num_invalid_calls_fuzzy = sv_db.get_num_invalid_calls(id_a, 0, blur_amount)
    avg_blur = sv_db.get_blur_on_overlaps_between_calls(id_b, id_a, 0, blur_amount)

    return xs_2, ys_2, ps, num_invalid_calls_fuzzy, avg_blur
    #return xs, ys, xs_2, ys_2, ps, num_invalid_calls, num_invalid_calls_fuzzy, avg_blur

def compare_caller(sv_db, id_a, id_b, min_score):
    num_calls_a = sv_db.get_num_calls(id_a, min_score) # num calls made
    num_calls_b = sv_db.get_num_calls(id_b, min_score) # num actual calls
    if num_calls_b == 0:
        return (0, 0, 0, 0, 0, 0)
    call_area_a = sv_db.get_call_area(id_a, min_score)
    if num_calls_a > 0:
        rel_call_area_a = int(math.sqrt(call_area_a/num_calls_a)) # get the edge length
    else:
        rel_call_area_a = "n/a"
    call_area_b = sv_db.get_call_area(id_b, min_score)
    rel_call_area_b = call_area_b/num_calls_b
    num_overlaps_a_to_b = sv_db.get_num_overlaps_between_calls(id_a, id_b, min_score, 0) # true positives
    num_almost_overlaps_a_to_b = sv_db.get_num_overlaps_between_calls(id_a, id_b, min_score, blur_amount) # true positives
    num_almost_overlaps_a_to_b -= num_overlaps_a_to_b
    num_overlaps_b_to_a = sv_db.get_num_overlaps_between_calls(id_b, id_a, min_score, 0) # how many of the sv's are detected?
    num_almost_overlaps_b_to_a = sv_db.get_num_overlaps_between_calls(id_b, id_a, min_score, blur_amount) # how many of the sv's are detected?
    num_errors = num_calls_b - num_overlaps_b_to_a # how many of the sv's are NOT detected?
    num_way_off = num_calls_b - num_almost_overlaps_b_to_a # how many of the sv's are NOT detected?
    num_invalid_calls = sv_db.get_num_invalid_calls(id_a, min_score, 0)
    return (num_calls_a, num_overlaps_b_to_a, num_almost_overlaps_b_to_a, num_errors, num_way_off, 
            rel_call_area_a, num_invalid_calls)

def compare_callers(db_name, names_a, names_b=["simulated sv"], min_scores=[0]):
    sv_db = SV_DB(db_name, "open")
    #print("sensitivity = true positive rate = recall")
    #print("missing rate = how many calls are missing")
    out = [["test caller", "ground truth caller", "min score", "#calls", "#found", "#almost", "#missed", "#way off", "fuzziness", "#inv call"]]
    for name_a, name_b in zip(names_a, names_b):
        id_a = sv_db.get_run_id(name_a)
        id_b = sv_db.get_run_id(name_b)
        date_a = sv_db.get_run_date(id_a) 
        date_b = sv_db.get_run_date(id_b)
        for min_score in min_scores:
            out.append([name_a + " - " + date_a, name_b + " - " + date_b, min_score,
                       *(str(x) for x in compare_caller(sv_db, id_a, id_b, min_score))])
    print_columns(out)

def compare_all_callers_against(sv_db, json_info_file, out_file_name=None, outfile_2_name=None):
    out = []
    out_2 = {}
    for dataset in json_info_file["datasets"]:
        id_b = dataset["ground_truth"]
        name_b = dataset["name"]
        date_b = sv_db.get_run_date(id_b)
        print("ground truth is:", name_b, "-", date_b, "[ id:", id_b, "] - with", sv_db.get_num_calls(id_b, 0), "calls")
        #print("sensitivity = true positive rate = recall")
        #print("missing rate = how many calls are missing")
        for id_a in sv_db.newest_unique_runs(2, "ground_truth=" + str(id_b)):
            name_a = sv_db.get_run_name(id_a).split("-")
            dataset_name = name_a[0]
            dataset_size = name_a[1]
            seq = name_a[2]
            cov = name_a[3]
            aligner = name_a[4]
            caller = name_a[5]
            date_a = sv_db.get_run_date(id_a)
            print("analyzing", id_a, name_a, date_a)
            out.append([dataset_name, dataset_size, seq, cov, caller, aligner, str(id_a),
                        *(str(x) for x in compare_caller(sv_db, id_a, id_b, 60 if caller == "MA_SV" else 0))])
            if str(name_a[:4]) not in out_2:
                out_2[str(name_a[:4])] = {}
            out_2[str(name_a[:4])][str((aligner, caller))] = analyze_by_score(sv_db, id_a, id_b)

    out.sort()
    out.insert(0, ["dataset", "size", "sequencer", "coverage", "caller", "aligner", "id", "#calls", "#found", "#almost",
                    "#missed", "#way off", "fuzziness", "#inv calls"])
    print()
    print_columns(out)
    if not out_file_name is None:
        with open(out_file_name, "w") as file_out:
            for line in out:
                for cell in line:
                    file_out.write(cell)
                    file_out.write("\t")
                file_out.write("\n")
    if not outfile_2_name is None:
        with open(outfile_2_name, "w") as file_out:
            json.dump(out_2, file_out)


def analyze_sample_dataset(dataset_name, run_callers=True, recompute_jumps=False, out_file_name=None):
    # decode hook for the json that decodes lists dicts and floats properly
    def _decode(o):
        if isinstance(o, str):
            try:
                return float(o)
            except ValueError:
                return o
        elif isinstance(o, dict):
            return {_decode(k): _decode(v) for k, v in o.items()}
        elif isinstance(o, list):
            return [_decode(v) for v in o]
        else:
            return o
    #actually open and load the info.json file
    json_info_file = None # noop
    with open("/MAdata/sv_datasets/" + dataset_name + "/info.json", "r") as json_file:
        json_info_file = json.loads(json_file.read(), object_hook=_decode)

    # create the calls
    db = SV_DB("/MAdata/sv_datasets/" + dataset_name + "/svs.db", "open")
    if run_callers:
        pack = Pack()
        pack.load(json_info_file["reference_path"] + "/ma/genome")
        fm_index = FMIndex()
        fm_index.load(json_info_file["reference_path"] + "/ma/genome")

        # create alignment files if they do not exist
        create_alignments_if_necessary(dataset_name, json_info_file, db, pack, fm_index, recompute_jumps)
        # save the info.json file
        print(json_info_file)
        with open("/MAdata/sv_datasets/" + dataset_name + "/info.json", "w") as json_out:
            json.dump(json_info_file, json_out)


        run_callers_if_necessary(dataset_name, json_info_file, db, pack, fm_index)

        # save the info.json file
        print(json_info_file)
        with open("/MAdata/sv_datasets/" + dataset_name + "/info.json", "w") as json_out:
            json.dump(json_info_file, json_out)

    compare_all_callers_against(db, json_info_file, "/MAdata/sv_datasets/" + dataset_name + "/bar_diagrams.tsv", "/MAdata/sv_datasets/" + dataset_name + "/by_score.json")


#print("===============")
#compare_callers("/MAdata/databases/sv_simulated", ["MA-SV"])
#print("===============")
if __name__ == "__main__":
    analyze_sample_dataset("minimal", True)
    #analyze_sample_dataset("comprehensive_random", False)
    
    #compare_all_callers_against(SV_DB("/MAdata/databases/sv_simulated", "open"))