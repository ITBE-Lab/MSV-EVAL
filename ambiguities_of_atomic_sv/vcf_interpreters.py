from MS import *
from MA import *
from MSV import *
import math
import traceback
import os

logged_errors = set()
def log_error(call, error_file, interpreter_name, e=None):
    # don't log an error twice
    key = interpreter_name + call["ALT"]
    if "SVTYPE" in call["INFO"]:
        key = interpreter_name + call["INFO"]["SVTYPE"]
    if key in logged_errors:
        return
    logged_errors.add(key)

    print("unrecognized sv:", call)
    error_file.write("============== unrecognized sv ==============\n")
    error_file.write("in interpreter: " + interpreter_name + "\n")
    error_file.write(str(call))
    if not e is None:
        error_file.write("\n")
        error_file.write(str(e))
        error_file.write(traceback.format_exc())
    error_file.write("\n\n\n")

def sniffles_interpreter(call, pack, error_file):
    def find_confidence(call):
        if call["FILTER"] != "PASS":
            return 0
        return int(float(call["INFO"]["RE"]))

    def find_from_and_to_pos(call):
        from_pos = int(call["POS"]) + pack.start_of_sequence(call["CHROM"])
        to_pos = int(call["INFO"]["END"]) + pack.start_of_sequence(call["INFO"]["CHR2"])
        return from_pos, to_pos

    def find_std_from_std_to(call):
        std_from = math.ceil(float(call["INFO"]["STD_quant_start"]))
        std_to = math.ceil(float(call["INFO"]["STD_quant_stop"]))
        return std_from, std_to

    try:
        from_pos, to_pos = find_from_and_to_pos(call)
        if "PRECISE" in call["INFO"]:
            std_from, std_to = (0, 0)
        elif "IMPRECISE" in call["INFO"]:
            std_from, std_to = find_std_from_std_to(call)
            #underflow protection
            if from_pos < std_from//2:
                std_from = from_pos//2
            #underflow protection
            if to_pos < std_to//2:
                std_to = to_pos//2
            from_pos -= std_from//2
            to_pos -= std_to//2
        else:
            raise Exception("found neither precise nor imprecise in INFO")

        return from_pos, to_pos, int(call["ID"]), call["ALT"] + "-conf:" + str(find_confidence(call))
    except Exception as e:
        log_error(call, error_file, "sniffles", e)

def delly_interpreter(call, pack, error_file):
    def find_confidence(call):
        if call["FILTER"] != "PASS":
            return 0
        conf = 0
        if "PE" in call["INFO"]:
            conf += int(call["INFO"]["PE"])
        if "SR" in call["INFO"]:
            conf += int(call["INFO"]["SR"])
        return conf

    def find_std_from_std_to(call):
        std_from = call["INFO"]["CIPOS"].split(",")
        std_to = call["INFO"]["CIEND"].split(",")
        return math.ceil(float(std_from[1]) - float(std_from[0])), math.ceil(float(std_to[1]) - float(std_to[0]))

    def find_from_and_to_pos(call):
        from_pos = int(call["POS"]) + pack.start_of_sequence(call["CHROM"]) - 1
        if "CHR2" in call["INFO"]:
            to_pos = int(call["INFO"]["END"]) + pack.start_of_sequence(call["INFO"]["CHR2"]) - 1
        else:
            to_pos = int(call["INFO"]["END"]) + pack.start_of_sequence(call["CHROM"]) - 1
        return from_pos, to_pos

    try:
        from_pos, to_pos = find_from_and_to_pos(call)
        if "PRECISE" in call["INFO"]:
            std_from, std_to = (0, 0)
        elif "IMPRECISE" in call["INFO"]:
            std_from, std_to = find_std_from_std_to(call)
            #underflow protection
            if from_pos < std_from//2:
                std_from = from_pos//2
            #underflow protection
            if to_pos < std_to//2:
                std_to = to_pos//2
            from_pos -= int(std_from/2)
            to_pos -= int(std_to/2)
        else:
            raise Exception("found neither precise nor imprecise in INFO")

        call_name = call["ALT"] + " " + call["INFO"]["SVTYPE"]
        if "CT" in call["INFO"]:
            call_name += " " + call["INFO"]["CT"]

        return from_pos, to_pos, int(call["ID"][4:]), call_name + "-conf:" + str(find_confidence(call))
    except Exception as e:
        log_error(call, error_file, "delly", e)

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
        for line in vcf_file:
            if line[-1] == "\n":
                line = line[:-1]
            if line[:2] == "##":
                info.append(line)
            elif line[0] == "#":
                names = line[1:].split("\t")
            else:
                d = {}
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