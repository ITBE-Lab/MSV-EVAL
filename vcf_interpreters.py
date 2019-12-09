from MA import *
import math
import traceback

bnd_mate_dict = {}
def default_vcf_interpreter(call, call_inserter, pack, error_file):
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

def sniffles_vcf_interpreter(call, call_inserter, pack, error_file):
    if call["ALT"] == "<DEL>" and "PRECISE" in call["INFO"]:
        #print(call)
        from_pos = int(call["POS"]) + pack.start_of_sequence(call["CHROM"])
        to_pos = int(call["INFO"]["END"]) + pack.start_of_sequence(call["INFO"]["CHR2"])
        call_inserter.insert_call(SvCall(from_pos, to_pos, 0, 0, False, int(float(call["INFO"]["RE"])), 1))
    elif call["ALT"] == "<DEL>" and "IMPRECISE" in call["INFO"]:
        #print(call)
        std_from, std_to = find_std_from_std_to(call)
        from_pos = int(call["POS"]) + pack.start_of_sequence(call["CHROM"]) - int(std_from/2)
        if "CHR2" in call["INFO"]:
            to_pos = int(call["INFO"]["END"]) + pack.start_of_sequence(call["INFO"]["CHR2"]) - int(std_to/2)
        else:
            to_pos = from_pos - int(call["INFO"]["SVLEN"]) - int(std_to/2)
        call_inserter.insert_call(SvCall(from_pos, to_pos, std_from, std_to, False, int(float(call["INFO"]["RE"])), 1))
    else:
        print("unrecognized sv:", call)

def pb_sv_interpreter(call, call_inserter, pack, error_file):
    def find_confidence(call):
        idx = call["FORMAT"].split(":").index("SAC")
        return sum([int(x) for x in call[call.names[-1]].split(":")[idx].split(",")])

    if call["INFO"]["SVTYPE"] == "DEL":
        from_pos = int(call["POS"]) + pack.start_of_sequence(call["CHROM"])
        to_pos = from_pos - int(call["INFO"]["SVLEN"])
        call_inserter.insert_call(SvCall(from_pos, to_pos, 0, 0, False, find_confidence(call), 1))
    else:
        print("unrecognized sv:", call)

def delly_interpreter(call, call_inserter, pack, error_file):
    def find_confidence(call):
        return int(call["INFO"]["PE"]) + int(call["INFO"]["SR"])

    def find_std_from_std_to(call):
        std_from = call["INFO"]["CIPOS"].split(",")
        std_to = call["INFO"]["CIEND"].split(",")
        return math.ceil(float(std_from[1]) - float(std_from[0])), math.ceil(float(std_to[1]) - float(std_to[0]))

    if call["ALT"] == "<DEL>" and "PRECISE" in call["INFO"]:
        #print(call)
        from_pos = int(call["POS"]) + pack.start_of_sequence(call["CHROM"])
        to_pos = int(call["INFO"]["END"]) + pack.start_of_sequence(call["INFO"]["CHR2"])
        call_inserter.insert_call(SvCall(from_pos, to_pos, 0, 0, False, find_confidence(call), 1))
    elif call["ALT"] == "<DEL>" and "IMPRECISE" in call["INFO"]:
        #print(call)
        std_from, std_to = find_std_from_std_to(call)
        from_pos = int(call["POS"]) + pack.start_of_sequence(call["CHROM"]) - int(std_from/2)
        to_pos = int(call["INFO"]["END"]) + pack.start_of_sequence(call["INFO"]["CHR2"]) - int(std_to/2)
        call_inserter.insert_call(SvCall(from_pos, to_pos, std_from, std_to, False, find_confidence(call), 1))
    else:
        print("unrecognized sv:", call)