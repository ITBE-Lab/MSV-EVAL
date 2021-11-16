from MS import *
from MA import *
from MSV import *
import math
import traceback
import os
import fnmatch
from sv_util.os_aligners import bwa, sam_to_bam
from sv_util.settings import *
from ambiguities_of_atomic_sv.os_sv_callers import gridss
from ambiguities_of_atomic_sv.vcf_interpreters import log_error

#genome_dir = global_prefix + "genome/human/GRCh38.p12/"
genome_dir = main_data_folder + "/genome/yeasts/YPS138/"
read_data_dir = main_data_folder + "/ena/" #"giab/"

def regex_match(folder, regex):
    l = []
    for file in os.listdir(folder):
        if fnmatch.fnmatch(file, regex):
            l.append(folder + str(file))
    return l

bnd_mate_dict_gridss = {}
def gridss_interpreter(call, pack, error_file):
    def find_confidence(call):
        if call["FILTER"] != "PASS":
            return 0
        return int(float(call["QUAL"]))

    def does_call_reverse_strand(call):
        if call["ALT"][-1] == "[":
            return False
        if call["ALT"][-1] == "]":
            return True
        elif call["ALT"][0] == "[":
            return True
        elif call["ALT"][0] == "]":
            return False
        else:
            raise Exception("could not classify call")
    def does_call_come_from_reverse_strand(call):
        if call["ALT"][-1] == "[":
            return True
        if call["ALT"][-1] == "]":
            return False
        elif call["ALT"][0] == "[":
            return True
        elif call["ALT"][0] == "]":
            return False
        else:
            raise Exception("could not classify call")

    def get_insertion(call):
        if call["ALT"][0] in "[":
            return call["ALT"][call["ALT"].rfind("["):]
        elif call["ALT"][0] in "]":
            return call["ALT"][call["ALT"].rfind("]"):]
        elif call["ALT"][-1] in "[":
            return call["ALT"][:call["ALT"].find("[")]
        elif call["ALT"][-1] in "]":
            return call["ALT"][:call["ALT"].find("]")]
        else:
            raise Exception("could not classify call")

    try:
        if call["FILTER"] != "PASS":
            return None, ""
        if call["INFO"]["SVTYPE"] == "BND":
            if "MATEID" in call["INFO"]:
                if call["INFO"]["MATEID"] in bnd_mate_dict_gridss:
                    ins = get_insertion(call)
                    mate = bnd_mate_dict_gridss[call["INFO"]["MATEID"]]
                    from_pos = int(mate["POS"]) + pack.start_of_sequence(mate["CHROM"])
                    to_pos = int(call["POS"]) + pack.start_of_sequence(call["CHROM"])
                    assert does_call_reverse_strand(call) == does_call_reverse_strand(mate)
                    assert does_call_come_from_reverse_strand(call) != does_call_come_from_reverse_strand(mate)
                    from_forw = not does_call_come_from_reverse_strand(call)
                    to_forw = not from_forw if does_call_reverse_strand(call) else from_forw
                    return SvJump(from_pos, to_pos, 0, len(ins), from_forw, to_forw, find_confidence(call), \
                                  call["line_idx"], -1), ins
                    del bnd_mate_dict_gridss[call["INFO"]["MATEID"]]
                else:
                    bnd_mate_dict_gridss[call["ID"]] = call
                    return None, ""
            else:
                from_pos = int(call["POS"]) + pack.start_of_sequence(call["CHROM"])
                ins = get_insertion(call)
                return SvJump(from_pos, from_pos, 0, len(ins), True, True, find_confidence(call), \
                              call["line_idx"], -1), ins
        else:
            raise Exception("could not classify call")

    except Exception as e:
        log_error(call, error_file, "gridss", e)

def vcf_parser(file_name):
    class VCFFile:
        def __init__(self, d, names, layer, info):
            self.data = d
            self.names = names
            self.layer = layer
            self.info = info

        def __getitem__(self, name):
            if not name in self.data:
                return []
            return self.data[name]

        def __contains__(self, name):
            return name in self.data

        def __str__(self):
            s = ""
            #for info_line in self.info:
            #    s += info_line + "\n"
            s += "{\n"
            for key, val in self.data.items():
                for _ in range(self.layer + 1):
                    s += "\t"
                s += str(key) + ": " + str(val) + "\n"
            for _ in range(self.layer):
                s += "\t"
            return s + "}"

        def from_format(self, key, value_list_idx=-1):
            idx = self["FORMAT"].split(":").index(key)
            return self[self.names[value_list_idx]].split(":")[idx]

        def has_format(self, key):
            return key in self["FORMAT"].split(":")

    with open(file_name, "r") as vcf_file:
        names = []
        info = []
        for idx, line in enumerate(vcf_file):
            if line[-1] == "\n":
                line = line[:-1]
            if line[:2] == "##":
                info.append(line)
            elif line[0] == "#":
                names = line[1:].split("\t")
            else:
                d = {}
                d["line_idx"] = idx
                for name, field in zip(names, line.split("\t")):
                    if name == "INFO":
                        d2 = {}
                        keys = []
                        for key_value in field.split(";"):
                            if "=" in key_value:
                                key, value = key_value.split("=")
                                keys.append(key)
                                d2[key] = value
                            else:
                                d2[key_value] = True
                        d[name] = VCFFile(d2, keys, 1, [])
                    else:
                        d[name] = field
                yield VCFFile(d, names, 0, info)

def parse_and_insert(file_name, db_name, reference_genome, caller_name="gridss"):
    parameter_set_manager = ParameterSetManager()
    parameter_set_manager.set_selected("SV-PacBio")
    min_size = parameter_set_manager.by_name("Min Size Edge").get()
    db_conn = DbConn({"SCHEMA": {"NAME": db_name}})

    
    pack = Pack()
    pack.load(reference_genome)

    JumpRunTable(db_conn)
    SvCallerRunTable(db_conn)
    sv_call_table = SvCallTable(db_conn)
    #call_per_contig_table = FirstCallPerContigTable(db_conn)
    caller_run_id = GetCallInserter(parameter_set_manager, db_conn, caller_name + " - Large",
                                   "SV's form genome assembly that are too long for Illumina reads", -1).cpp_module.id
    contig_border_caller_run_id = GetCallInserter(parameter_set_manager, db_conn, caller_name + " - Contig Border",
                                   "SV's form genome assembly that are not covered by any alignments", -1).cpp_module.id
    small_caller_run_id = GetCallInserter(parameter_set_manager, db_conn, caller_name + " - Small",
                                   "SV's form genome assembly that are too small for pacBio reads", -1).cpp_module.id

    contig_filter = JumpsFilterContigBorder(parameter_set_manager)

    cnt_small = 0
    cnt_large = 0
    cnt_contig_border = 0
    with open(file_name + ".vcf_errors.log", "w") as error_file:
        for vcf_call in vcf_parser(file_name):
            jump, ins = gridss_interpreter(vcf_call, pack, error_file)
            if not jump is None:
                f = jump.from_pos
                t = jump.to_pos
                if (not jump.from_known() or not jump.to_known()) and not jump.from_forward:
                    # flip from and to information for dummy jumps from reverse strand seeds
                    # Why do we need to do this?:
                    # On reads we do not know the orientation of the reads. So, in order to unify forward and reverse
                    # strand reads, we mirror the dummy calls as well. Since they technically do not have a x and y position
                    # they cannot get the mirrored flag. Instead we simply mirror all reverse strand dummy calls.
                    # Hence these calls are wrong in the context of the reconstruction.
                    # since we do not (and cannot since there is no coverage at chrom ends) analyze these calls we can
                    # place them as we wish... we choose to simply mirror everything back to the original position
                    # since i do not want to do this in the cpp code i do it here in the py code.
                    tmp = f
                    f = t
                    t = tmp

                call = SvCall(f, t, 0, 0, jump.from_forward, jump.to_forward, 1, jump.num_supp_nt())
                call.inserted_sequence = NucSeq(ins)
                call.inserted_sequence.check()
                call.id = jump.id
                call.mirrored = jump.was_mirrored
                if contig_filter.cpp_module.by_contig_border(jump, pack):
                    cnt_contig_border += 1
                    sv_call_table.insert_call(contig_border_caller_run_id, call)
                elif max(abs(f - t), len(ins)) < min_size:
                    cnt_small += 1
                    sv_call_table.insert_call(small_caller_run_id, call)
                else:
                    cnt_large += 1
                    sv_call_table.insert_call(caller_run_id, call)


    print("Inserted into DB. There are", cnt_small, "small and", cnt_large, "large entries.", cnt_contig_border,
          "entries are too close to contig borders and are filtered out.")

    return caller_run_id

db_name = "UFRJ50816"

def main():
    reads = regex_match(read_data_dir + "simulated/UFRJ50816/Illumina-250/", "*.bwa.read1.fastq.gz")
    reads_mates = regex_match(read_data_dir + "simulated/UFRJ50816/Illumina-250/", "*.bwa.read2.fastq.gz")

    json_dict = {"reference_path":genome_dir}
    read_json = {"technology":"illumina", "name":"n/a", "fasta_file":", ".join(reads),
                 "fasta_file_mate":", ".join(reads_mates)}
    path_sam = gridss_data_dir + ".bwa_alignment"
    #bwa(read_json, path_sam + ".sam", json_dict)
    #sam_to_bam(path_sam)
    #gridss(path_sam + ".sorted.bam", gridss_data_dir+".gridss.vcf", genome_dir+"fasta/genome.fna")

    run_id = parse_and_insert(gridss_data_dir+".gridss.vcf", db_name, genome_dir+"ma/genome")
    print("run_id", run_id)

main()