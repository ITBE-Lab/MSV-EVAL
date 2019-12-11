from MA import *
import math
import traceback

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

bnd_mate_dict = {}
def default_vcf_interpreter(call, call_inserter, pack, error_file):
    def find_confidence(call):
        if call["FILTER"] != "PASS":
            return 0
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
        return 1
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
    try:
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
            if "CHR2" in call["INFO"]:
                to_pos = int(call["INFO"]["END"]) + pack.start_of_sequence(call["INFO"]["CHR2"]) - int(std_to/2)
            else:
                to_pos = from_pos - int(call["INFO"]["SVLEN"]) - int(std_to/2)
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
            to_pos = int(call["INFO"]["END"]) + pack.start_of_sequence(call["CHROM"])
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
                call_inserter.insert_call(SvCall(from_pos, to_pos, std_to, std_from, save_as_inversion,
                                                    find_confidence(call), 1))
                del bnd_mate_dict[call["INFO"]["MATEID"]]
            else:
                bnd_mate_dict[call["ID"]] = call
        else:
            log_error(call, error_file, "default")
            #exit(0)
    except Exception as e:
        log_error(call, error_file, "default", e)

def sniffles_interpreter(call, call_inserter, pack, error_file):
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

        to_recognize = 1
        if "/" in call["ALT"]:
            to_recognize = 2
        if "DEL" in call["ALT"]:
            call_inserter.insert_call(SvCall(from_pos, to_pos, std_from, std_to, False, find_confidence(call), 1))
            to_recognize -= 1
        if "INV" in call["ALT"]:
            call_inserter.insert_call(SvCall(from_pos, to_pos, std_from, std_to, True, find_confidence(call), 1))
            call_inserter.insert_call(SvCall(to_pos, from_pos, std_from, std_to, True, find_confidence(call), 1))
            to_recognize -= 1
        if "INS" in call["ALT"]:
            call_inserter.insert_call(SvCall(from_pos, from_pos + 1, std_from, std_from, False,
                                      find_confidence(call), 1))
            to_recognize -= 1
        if "DUP" in call["ALT"]:
            call_inserter.insert_call(SvCall(to_pos, from_pos, std_to, std_from, False, find_confidence(call), 1))
            to_recognize -= 1
        if "TRA" in call["ALT"]:
            call_inserter.insert_call(SvCall(from_pos, to_pos, std_from, std_to, False, find_confidence(call), 1))
            to_recognize -= 1

        if to_recognize != 0:
            raise Exception("could not classify call")
    except Exception as e:
        log_error(call, error_file, "sniffles", e)

bnd_mate_dict_pb_sv = {}
def pb_sv_interpreter(call, call_inserter, pack, error_file):
    def find_confidence(call):
        if call["FILTER"] != "PASS":
            return 0
        if call.has_format("SAC"):
            sac = call.from_format("SAC")
            if "," in sac:
                return sum([int(x) for x in sac.split(",")])
            return int(sac)
        elif call.has_format("DP"):
            return int(call.from_format("DP"))
        return 0
    def find_from_and_to_pos(call):
        from_pos = int(call["POS"]) + pack.start_of_sequence(call["CHROM"])
        if "END" in call["INFO"]:
            to_pos = int(call["INFO"]["END"]) + pack.start_of_sequence(call["CHROM"])
        else:
            to_pos = None
        return from_pos, to_pos
    def find_std_from(call):
        std_from = call["INFO"]["CIPOS"].split(",")
        return math.ceil(float(std_from[1]) - float(std_from[0]))

    try:
        from_pos, to_pos = find_from_and_to_pos(call)
        if call["INFO"]["SVTYPE"] == "DEL":
            call_inserter.insert_call(SvCall(from_pos, to_pos, 0, 0, False, find_confidence(call), 1))
            return
        if call["INFO"]["SVTYPE"] == "INS":
            call_inserter.insert_call(SvCall(from_pos, from_pos + 1, 0, 0, False, find_confidence(call), 1))
            return
        if call["INFO"]["SVTYPE"] == "INV":
            call_inserter.insert_call(SvCall(from_pos, to_pos, 0, 0, True, find_confidence(call), 1))
            call_inserter.insert_call(SvCall(to_pos, from_pos, 0, 0, True, find_confidence(call), 1))
            return
        if call["INFO"]["SVTYPE"] == "DUP":
            call_inserter.insert_call(SvCall(to_pos, from_pos, 0, 0, False, find_confidence(call), 1))
            return
        if call["INFO"]["SVTYPE"] == "cnv":
            return
        if call["INFO"]["SVTYPE"] == "BND":
            if call["INFO"]["MATEID"] in bnd_mate_dict_pb_sv:
                mate = bnd_mate_dict_pb_sv[call["INFO"]["MATEID"]]
                std_from = find_std_from(call)
                std_to = find_std_from(mate)
                from_pos = int(mate["POS"]) + pack.start_of_sequence(mate["CHROM"]) - std_from//2
                to_pos = int(call["POS"]) + pack.start_of_sequence(call["CHROM"]) - std_to//2
                call_inserter.insert_call(SvCall(from_pos, to_pos, from_pos, std_to, False, find_confidence(call), 1))
                call_inserter.insert_call(SvCall(to_pos, from_pos, std_to, from_pos, False, find_confidence(mate), 1))
                del bnd_mate_dict_pb_sv[call["INFO"]["MATEID"]]
            else:
                bnd_mate_dict_pb_sv[call["ID"]] = call
            return
        raise Exception("could not classify call")
    except Exception as e:
        log_error(call, error_file, "pb_sv", e)

def delly_interpreter(call, call_inserter, pack, error_file):
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
        to_pos = int(call["INFO"]["END"]) + pack.start_of_sequence(call["INFO"]["CHR2"]) - 1
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

        if call["ALT"] == "<DEL>":
            call_inserter.insert_call(SvCall(from_pos, to_pos, std_from, std_to, False, find_confidence(call), 1))
            return
        if call["ALT"] == "<INV>" and False:
            from_pos, to_pos = find_from_and_to_pos(call)
            # delly calls inversion twice once 3t3 once 5t5 (or heat to head and tail to tail) of reads
            call_inserter.insert_call(SvCall(from_pos, to_pos, 0, 0, True, find_confidence(call), 1))
            call_inserter.insert_call(SvCall(to_pos, from_pos, 0, 0, True, find_confidence(call), 1))
            return
        if call["ALT"] == "<INV>":
            # delly calls inversion twice once 3t3 once 5t5 (or heat to head and tail to tail) of reads
            call_inserter.insert_call(SvCall(from_pos, to_pos, std_from, std_to, True, find_confidence(call), 1))
            call_inserter.insert_call(SvCall(to_pos, from_pos, std_from, std_to, True, find_confidence(call), 1))
            return
        if call["ALT"] == "<DUP>":
            call_inserter.insert_call(SvCall(to_pos, from_pos, std_from, std_to, False, find_confidence(call), 1))
            return
        if call["ALT"] == "<INS>":
            call_inserter.insert_call(SvCall(from_pos, from_pos + 1, std_from, std_to, False, find_confidence(call), 1))
            return
        if call["INFO"]["SVTYPE"] == "BND":
            call_inserter.insert_call(SvCall(to_pos, from_pos, std_from, std_to, False, find_confidence(call), 1))
            call_inserter.insert_call(SvCall(from_pos, to_pos, std_to, std_from, False, find_confidence(call), 1))
            return
        raise Exception("could not classify call")

    except Exception as e:
        log_error(call, error_file, "delly", e)

bnd_mate_dict_manta = {}
def manta_interpreter(call, call_inserter, pack, error_file):
    def find_confidence(call):
        if call["FILTER"] != "PASS":
            return 0
        return int(call["QUAL"])

    def find_std_from_std_to(call):
        std_from = call["INFO"]["CIPOS"].split(",")
        if "CIEND" in call["INFO"]:
            std_to = call["INFO"]["CIEND"].split(",")
        else:
            std_to = (0, 0)
        return math.ceil(float(std_from[1]) - float(std_from[0])), math.ceil(float(std_to[1]) - float(std_to[0]))

    def find_from_and_to_pos(call):
        from_pos = int(call["POS"]) + pack.start_of_sequence(call["CHROM"])
        if "END" in call["INFO"]:
            to_pos = int(call["INFO"]["END"]) + pack.start_of_sequence(call["CHROM"])
        else:
            to_pos = 0
        return from_pos, to_pos

    try:
        from_pos, to_pos = find_from_and_to_pos(call)
        if "IMPRECISE" in call["INFO"]:
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
            std_from, std_to = (0, 0)

        if call["ALT"] == "<DUP:TANDEM>":
            call_inserter.insert_call(SvCall(to_pos, from_pos, std_from, std_to, False, find_confidence(call), 1))
            return
        if call["ALT"] == "<DUP>":
            call_inserter.insert_call(SvCall(to_pos, from_pos, std_from, std_to, False, find_confidence(call), 1))
            return
        if call["INFO"]["SVTYPE"] == "DEL":
            call_inserter.insert_call(SvCall(from_pos, to_pos, std_from, std_to, False, find_confidence(call), 1))
            return
        if call["INFO"]["SVTYPE"] == "BND":
            if call["INFO"]["MATEID"] in bnd_mate_dict_manta:
                mate = bnd_mate_dict_manta[call["INFO"]["MATEID"]]
                std_from = find_std_from(call)
                std_to = find_std_from(mate)
                from_pos = int(mate["POS"]) + pack.start_of_sequence(mate["CHROM"]) - std_from//2
                to_pos = int(call["POS"]) + pack.start_of_sequence(call["CHROM"]) - std_to//2
                call_inserter.insert_call(SvCall(from_pos, to_pos, from_pos, std_to, False, find_confidence(call), 1))
                call_inserter.insert_call(SvCall(to_pos, from_pos, std_to, from_pos, False, find_confidence(mate), 1))
                del bnd_mate_dict_manta[call["INFO"]["MATEID"]]
            else:
                bnd_mate_dict_manta[call["ID"]] = call
        raise Exception("could not classify call")

    except Exception as e:
        log_error(call, error_file, "manta", e)

def smoove_interpreter(call, call_inserter, pack, error_file):
    def find_confidence(call):
        return int(float(call["QUAL"])*100)

    def find_std_from_std_to(call):
        std_from = call["INFO"]["CIPOS"].split(",")
        if "CIEND" in call["INFO"]:
            std_to = call["INFO"]["CIEND"].split(",")
        else:
            std_to = (0, 0)
        return math.ceil(float(std_from[1]) - float(std_from[0])), math.ceil(float(std_to[1]) - float(std_to[0]))

    def find_from_and_to_pos(call):
        from_pos = int(call["POS"]) + pack.start_of_sequence(call["CHROM"])
        if "END" in call["INFO"]:
            to_pos = int(call["INFO"]["END"]) + pack.start_of_sequence(call["CHROM"])
        else:
            to_pos = 0
        return from_pos, to_pos

    try:
        from_pos, to_pos = find_from_and_to_pos(call)
        if "IMPRECISE" in call["INFO"]:
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
            std_from, std_to = (0, 0)

        if call["ALT"] == "<DUP:TANDEM>":
            call_inserter.insert_call(SvCall(to_pos, from_pos, std_from, std_to, False, find_confidence(call), 1))
            return
        if call["ALT"] == "<DUP>":
            call_inserter.insert_call(SvCall(to_pos, from_pos, std_from, std_to, False, find_confidence(call), 1))
            return
        if call["INFO"]["SVTYPE"] == "DEL":
            call_inserter.insert_call(SvCall(from_pos, to_pos, std_from, std_to, False, find_confidence(call), 1))
            return
        if call["INFO"]["SVTYPE"] == "BND":
            if call["INFO"]["MATEID"] in bnd_mate_dict_manta:
                mate = bnd_mate_dict_manta[call["INFO"]["MATEID"]]
                std_from = find_std_from(call)
                std_to = find_std_from(mate)
                from_pos = int(mate["POS"]) + pack.start_of_sequence(mate["CHROM"]) - std_from//2
                to_pos = int(call["POS"]) + pack.start_of_sequence(call["CHROM"]) - std_to//2
                call_inserter.insert_call(SvCall(from_pos, to_pos, from_pos, std_to, False, find_confidence(call), 1))
                call_inserter.insert_call(SvCall(to_pos, from_pos, std_to, from_pos, False, find_confidence(mate), 1))
                del bnd_mate_dict_manta[call["INFO"]["MATEID"]]
            else:
                bnd_mate_dict_manta[call["ID"]] = call
        raise Exception("could not classify call")

    except Exception as e:
        log_error(call, error_file, "smoove", e)