from MA import *
import math
import os
import json
import compute_sv_jumps
import sweep_sv_jumps

def create_alignments_if_necessary(dataset_name, json_dict, db, pack, fm_index):
    def bwa_paired(read_set, sam_file_path):
        index_str = json_dict["reference_path"] + "/bwa/genome"
        os.system("~/workspace/bwa/bwa mem -t 32 " + index_str + " " + read_set["fasta_file"] + " "
                  + read_set["fasta_file_mate"] + " > " + sam_file_path + " 2> /dev/null")

    def minimap2(read_set, sam_file_path):
        presetting = None # noop
        if read_set["technology"] == "pb":
            presetting = "map-pb"
        if read_set["technology"] == "ont":
            presetting = "map-ont"
        index_str = json_dict["reference_path"] + "/minimap/genome." + presetting + ".mmi"
        # -c output CIGAR in PAF; -a output SAM
        os.system("~/workspace/minimap2/minimap2 -c -a -t 32 -x " + presetting + " " + index_str + " "
                  + read_set["fasta_file"] + " > " + sam_file_path + " 2> /dev/null")

    def ngmlr(read_set, sam_file_path):
        presetting = None # noop
        if read_set["technology"] == "pb":
            presetting = "pacbio"
        if read_set["technology"] == "ont":
            presetting = "ont"
        index_str = json_dict["reference_path"] + "/ngmlr/genome.fna"
        # -c output CIGAR in PAF; -a output SAM
        os.system("~/workspace/ngmlr/bin/ngmlr-0.2.8/ngmlr -r " + index_str + " -q "
                  + read_set["fasta_file"] + " -t 32 -x " + presetting + " > " + sam_file_path + " 2> /dev/null")

    alignment_calls = {
        "create_reads_survivor": [minimap2, ngmlr],
        "create_illumina_reads_dwgsim": [bwa_paired]
    }

    for read_set in json_dict["create_reads_funcs"]:
        if not "alignments" in read_set:
            read_set["alignments"] = []
        if not "jump_id" in read_set:
            print("computing jumps for MA-SV on", read_set["name"])
            db.drop_jump_indices() # @todo should work with partial indices here

            read_set["jump_id"] = compute_sv_jumps.compute_sv_jumps(ParameterSetManager(), fm_index, pack, db, 
                                                                    read_set["seq_id"])
        for alignment_call in alignment_calls[read_set["func_name"]]:
            sam_file_path = "/MAdata/sv_datasets/" + dataset_name + "/alignments/" \
                          + read_set["name"] + "." + alignment_call.__name__
            if os.path.exists( sam_file_path + ".sam" ):
                continue
            if not alignment_call.__name__  in read_set["alignments"]:
                read_set["alignments"].append(alignment_call.__name__)
            else:
                print("recreating alignments!")
            print("creating alignment files for", read_set["name"], alignment_call.__name__)
            alignment_call(read_set, sam_file_path + ".sam")

            # create sorted and indexed bam files
            sam_tools_pref = "~/workspace/samtools/samtools "
            to_bam_cmd = sam_tools_pref + "view -Sb " + sam_file_path + ".sam > " + sam_file_path + ".bam"
            os.system(to_bam_cmd)
            sort_cmd = sam_tools_pref + "sort -@ 32 -m 1G " + sam_file_path + ".bam > " + sam_file_path + ".sorted.bam"
            os.system(sort_cmd + " 2> /dev/null")
            index_cmd = sam_tools_pref + "index " + sam_file_path + ".sorted.bam > " + sam_file_path + ".sorted.bam.bai"
            os.system(index_cmd)
    db.create_jump_indices() # @todo should work with partial indices here

def vcf_parser(file_name):
    with open(file_name, "r") as vcf_file:
        names = []
        for line in vcf_file:
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
                        d[name] = d2
                    else:
                        d[name] = field
                yield d

def vcf_to_db(name, sv_db, file_name, pack):
    sv_db.clear_calls_table_for_caller(name)
    call_inserter = SvCallInserter(sv_db, name, "no desc", -1) # -1 since there are no related sv jumps...
    num_calls = 0
    for call in vcf_parser(file_name):
        num_calls += 1
        if call["ALT"] == "<DEL>" and "PRECISE" in call["INFO"]:
            #print(call)
            from_pos = int(call["POS"]) + pack.start_of_sequence(call["CHROM"])
            to_pos = int(call["INFO"]["END"]) + pack.start_of_sequence(call["INFO"]["CHR2"])
            call_inserter.insert_call(SvCall(from_pos, to_pos, 1, 1, False, float('inf')))
        elif call["ALT"] == "<DEL>" and "IMPRECISE" in call["INFO"]:
            #print(call)
            std_from = math.ceil(float(call["INFO"]["STD_quant_start"]))
            std_to = math.ceil(float(call["INFO"]["STD_quant_stop"]))
            from_pos = int(call["POS"]) + pack.start_of_sequence(call["CHROM"]) - int(std_from/2)
            to_pos = int(call["INFO"]["END"]) + pack.start_of_sequence(call["INFO"]["CHR2"]) - int(std_to/2)
            call_inserter.insert_call(SvCall(from_pos, to_pos, std_from, std_to, False, float('inf')))
        elif call["ALT"] == "<DUP>" and "PRECISE" in call["INFO"]:
            #print(call)
            from_pos = int(call["POS"]) + pack.start_of_sequence(call["CHROM"])
            to_pos = int(call["INFO"]["END"]) + pack.start_of_sequence(call["INFO"]["CHR2"])
            call_inserter.insert_call(SvCall(to_pos, from_pos, 1, 1, False, float('inf')))
        elif call["ALT"] == "<INS>" and "PRECISE" in call["INFO"]:
            #print(call)
            from_pos = int(call["POS"]) + pack.start_of_sequence(call["CHROM"])
            call_inserter.insert_call(SvCall(from_pos, from_pos, 1, 1, False, float('inf')))
        elif call["ALT"] == "<INS>" and "IMPRECISE" in call["INFO"]:
            #print(call)
            from_start = int(call["POS"]) + pack.start_of_sequence(call["CHROM"])
            from_end = int(call["INFO"]["END"]) + pack.start_of_sequence(call["INFO"]["CHR2"])
            call_inserter.insert_call(SvCall(from_start, from_start, from_end - from_start,
                                                from_end - from_start, False, float('inf')))
        elif call["ALT"] == "<INV>" and "PRECISE" in call["INFO"]:
            #print(call)
            from_pos = int(call["POS"]) + pack.start_of_sequence(call["CHROM"])
            to_pos = int(call["INFO"]["END"]) + pack.start_of_sequence(call["INFO"]["CHR2"])
            call_inserter.insert_call(SvCall(from_pos, to_pos, 1, 1, True, float('inf')))
            call_inserter.insert_call(SvCall(to_pos, from_pos, 1, 1, True, float('inf')))
        else:
            print("unrecognized sv:", call)
            exit(0)
    print("number of calls:", num_calls)

def run_callers_if_necessary(dataset_name, json_dict, db, pack):
    def sniffles(bam_file, vcf_file):
        os.system("~/workspace/Sniffles/bin/sniffles-core-1.0.8/sniffles -t 32 -m " + bam_file + " -v " + vcf_file
                  + " >/dev/null 2>&1")
    sv_calls = [sniffles]
    for read_set in json_dict["create_reads_funcs"]:
        if not "calls" in read_set:
            read_set["calls"] = []
        for sv_call in sv_calls:
            if not sv_call.__name__ in read_set["calls"]:
                read_set["calls"].append(sv_call.__name__)
            
            # MA-SV
            if not "MA_SV" in read_set["calls"]:
                read_set["calls"].append("MA_SV")
            print("creating calls for", read_set["name"], "MA_SV")
            sweep_sv_jumps.sweep_sv_jumps(ParameterSetManager(), db, read_set["jump_id"], 
                                          pack.unpacked_size_single_strand, read_set["name"] + "." + "MA_SV")

            for alignment in read_set["alignments"]:
                vcf_file_path = "/MAdata/sv_datasets/" + dataset_name + "/calls/" \
                        + read_set["name"] + "." + alignment + "." + sv_call.__name__ + ".vcf"
                bam_file_path = "/MAdata/sv_datasets/" + dataset_name + "/alignments/" \
                        + read_set["name"] + "." + alignment + ".sorted.bam"
                if os.path.exists( vcf_file_path ):
                    print("not creating calls for", read_set["name"], alignment, sv_call.__name__)
                    continue
                print("creating calls for", read_set["name"], alignment, sv_call.__name__)
                sv_call(bam_file_path, vcf_file_path)
                vcf_to_db(read_set["name"] + "." + alignment + "." + sv_call.__name__, db, vcf_file_path, pack)


def print_columns(data):
    col_width = [max([len(data[j][i]) for j in range(len(data))]) for i in range(len(data[0]))]
    first = True
    for row in data:
        print("| " + "".join(word.ljust(col_width[i]) + " | " for i, word in enumerate(row)))
        if first:
            print("-" * (sum(col_width) + len(col_width)*3 + 1))
            first = False

def compare_caller(sv_db, id_a, id_b, min_score):
    num_calls_a = sv_db.get_num_calls(id_a, min_score) # num calls made
    num_calls_b = sv_db.get_num_calls(id_b, min_score) # num actual calls
    call_area_a = sv_db.get_call_area(id_a, min_score)
    if num_calls_a > 0:
        rel_call_area_a = math.sqrt(call_area_a/num_calls_a) # get the edge length
    else:
        rel_call_area_a = 0
    call_area_b = sv_db.get_call_area(id_b, min_score)
    rel_call_area_b = call_area_b/num_calls_b
    num_overlaps_a_to_b = sv_db.get_num_overlaps_between_calls(id_a, id_b, min_score, 0) # true positives
    num_almost_overlaps_a_to_b = sv_db.get_num_overlaps_between_calls(id_a, id_b, min_score, 100) # true positives
    num_almost_overlaps_a_to_b -= num_overlaps_a_to_b
    num_overlaps_b_to_a = sv_db.get_num_overlaps_between_calls(id_b, id_a, min_score, 0) # how many of the sv's are detected?
    num_almost_overlaps_b_to_a = sv_db.get_num_overlaps_between_calls(id_b, id_a, min_score, 100) # how many of the sv's are detected?
    num_errors = num_calls_b - num_overlaps_b_to_a # how many of the sv's are NOT detected?
    return (num_calls_a, num_overlaps_b_to_a, num_almost_overlaps_b_to_a, num_errors, rel_call_area_a)

def compare_callers(db_name, names_a, names_b=["simulated sv"], min_scores=[0]):
    sv_db = SV_DB(db_name, "open")
    #print("sensitivity = true positive rate = recall")
    #print("missing rate = how many calls are missing")
    out = [["test caller", "ground truth caller", "min score", "#calls", "#found", "#almost", "#missed", "fuzziness"]]
    for name_a, name_b in zip(names_a, names_b):
        id_a = sv_db.get_run_id(name_a)
        id_b = sv_db.get_run_id(name_b)
        date_a = sv_db.get_run_date(id_a)
        date_b = sv_db.get_run_date(id_b)
        for min_score in min_scores:
            out.append([name_a + " - " + date_a, name_b + " - " + date_b, min_score,
                       *(str(x) for x in compare_caller(sv_db, id_a, id_b, min_score))])
    print_columns(out)

def compare_all_callers_against(sv_db, name_b="simulated sv"):
    id_b = sv_db.get_run_id(name_b)
    date_b = sv_db.get_run_date(id_b)
    #print("sensitivity = true positive rate = recall")
    #print("missing rate = how many calls are missing")
    print("ground truth is:", name_b, "-", date_b, "[ id:", id_b, "]")
    out = [["id", "test set", "time", "#calls", "#found", "#almost", "#missed", "fuzziness"]]
    for id_a in sv_db.newest_unique_runs(3):
        if id_a == id_b:
            continue
        name_a = sv_db.get_run_name(id_a)
        date_a = sv_db.get_run_date(id_a)
        out.append([str(id_a), name_a, date_a, *(str(x) for x in compare_caller(sv_db, id_a, id_b, 0))])
    print_columns(out)


def analyze_sample_dataset(dataset_name, run_callers=False):
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
        create_alignments_if_necessary(dataset_name, json_info_file, db, pack, fm_index)
        # save the info.json file
        print(json_info_file)
        with open("/MAdata/sv_datasets/" + dataset_name + "/info.json", "w") as json_out:
            json.dump(json_info_file, json_out)



        run_callers_if_necessary(dataset_name, json_info_file, db, pack)

        # save the info.json file
        print(json_info_file)
        with open("/MAdata/sv_datasets/" + dataset_name + "/info.json", "w") as json_out:
            json.dump(json_info_file, json_out)

    compare_all_callers_against(db)


#print("===============")
#compare_callers("/MAdata/databases/sv_simulated", ["MA-SV"])
#print("===============")
if __name__ == "__main__":
    analyze_sample_dataset("small_test_1", False)
    
    #compare_all_callers_against(SV_DB("/MAdata/databases/sv_simulated", "open"))