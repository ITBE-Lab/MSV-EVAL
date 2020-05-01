from MS import *
from MA import *
from MSV import *
import math
import os
from ..os_aligners import *
import json
import traceback
from pathlib import Path
import datetime
from bokeh.plotting import figure, show
from bokeh.models.formatters import PrintfTickFormatter
import caller_analysis.vcf_interpreters
from caller_analysis.os_sv_callers import *

# Setup the appropriate environment
"""Markus @ Zeus""" 
global_prefix = "/MAdata/"
svdb_dir = global_prefix + "sv_datasets3/"
sv_data_dir = global_prefix + "sv_datasets/"

"""Arne @ home"""
#global_prefix = "C:/MAdata/"
# svdb_dir = global_prefix + "sv_datasets/" 
# svdb_dir = global_prefix + "sv_datasets2/"

def create_alignments_if_necessary(dataset_name, json_dict, pack, fm_index, recompute_jumps=False, run_ma=True,
                                   run_others=True):
    alignment_calls = {
        "create_reads_survivor": [mm2, ngmlr, pbmm2],
        "create_illumina_reads_dwgsim": [bwa, bowtie]
    }

    with open(sv_data_dir + dataset_name + "/runtimes.log", "a") as runtime_file:
        for dataset in json_dict["datasets"]:
            for read_set in dataset["create_reads_funcs"]:
                if not "alignments" in read_set:
                    read_set["alignments"] = []
                if not "jump_id" in read_set or recompute_jumps:
                    if run_ma:
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
                        runtime_file.write(str(datetime.datetime.now()) + " " + dataset_name + " ")
                        runtime_file.write(dataset["name"] + " " + read_set["name"] + " compute_sv_jumps")
                        runtime_file.write("\n")
                        read_set["jump_id"] = compute_sv_jumps(params, fm_index, pack, dataset_name, read_set["seq_id"],
                                                                runtime_file)
                for alignment_call in alignment_calls[read_set["func_name"]]:
                    sam_file_path = sv_data_dir + dataset_name + "/alignments/" \
                                + read_set["name"] + "-" + alignment_call.__name__
                    if not alignment_call.__name__ in read_set["alignments"]:
                        read_set["alignments"].append(alignment_call.__name__)
                    if os.path.exists( sam_file_path + ".sam" ) or not run_others:
                        continue
                    print("creating alignment files for", read_set["name"], alignment_call.__name__)
                    alignment_call(read_set, sam_file_path + ".sam", json_dict)

                    if os.path.exists( sam_file_path + ".sam" ):
                        sam_to_bam(sam_file_path)
                    else:
                        print("aligner did not create alignments:", alignment_call)
def run_callers_if_necessary(dataset_name, json_dict, pack, fm_index, run_others=True, recompute_calls=False):
    # @todo svim?
    sv_calls = {
        "bwa":    [delly, manta], # smoove
        "bowtie": [delly, manta], # smoove
        "mm2":    [sniffles],
        "pbmm2":  [pbSv],
        "ngmlr":  [sniffles],
        "blasr":  [] #pbHoney @note the vcf file reading seems to be broken anyways
    }

    call_interpreters = {
        "delly": vcf_interpreters.delly_interpreter,
        "smoove": vcf_interpreters.smoove_interpreter,
        "manta": vcf_interpreters.manta_interpreter,
        "sniffles": vcf_interpreters.sniffles_interpreter,
        "pbSv": vcf_interpreters.pb_sv_interpreter
    }

    with open(sv_data_dir + dataset_name + "/vcf_errors.log", "a") as error_file:
        with open(sv_data_dir + dataset_name + "/runtimes.log", "a") as runtime_file:
            for dataset in json_dict["datasets"]:
                for read_set in dataset["create_reads_funcs"]:
                    if not "calls" in read_set:
                        read_set["calls"] = []

                    # MA-SV
                    if not "MA_SV" in read_set["calls"] or recompute_calls:
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
                        runtime_file.write(str(datetime.datetime.now()) + " " + dataset_name + " ")
                        runtime_file.write(dataset["name"] + " " + read_set["name"] + " sweep_sv_jumps_cpp")
                        runtime_file.write("\n")
                        if "jump_id" in read_set:
                            sweep_sv_jumps(params, dataset_name, read_set["jump_id"], read_set["name"] + "--" + "MA_SV",
                                            "ground_truth=" + str(dataset["ground_truth"]),
                                            [read_set["seq_id"]], pack, runtime_file)
                    # other callers
                    for alignment in read_set["alignments"]:
                        for sv_call in sv_calls[alignment]:
                            vcf_file_path = sv_data_dir + dataset_name + "/calls/" \
                                    + read_set["name"] + "-" + alignment + "-" + sv_call.__name__ + ".vcf"
                            bam_file_path = sv_data_dir + dataset_name + "/alignments/" \
                                    + read_set["name"] + "-" + alignment + ".sorted.bam"
                            if ( os.path.exists( vcf_file_path ) and not recompute_calls ) or not run_others:
                                print("not creating calls for", read_set["name"], alignment, sv_call.__name__)
                                continue
                            print("creating calls for", read_set["name"], alignment, sv_call.__name__)
                            sv_call(bam_file_path, vcf_file_path, json_dict["reference_path"])
                            if not os.path.exists( vcf_file_path ):
                                print("caller did not create calls: ", read_set["name"], alignment, sv_call.__name__)
                                continue
                            vcf_to_db(read_set["name"] + "-" + alignment + "-" + sv_call.__name__,
                                    "ground_truth=" + str(dataset["ground_truth"]), dataset_name, vcf_file_path, pack, 
                                    call_interpreters[sv_call.__name__], error_file,
                                    dataset["sv_func"] if "sv_func" in dataset else None)

                        for sv_call in sv_calls[alignment]:
                            if not sv_call.__name__ in read_set["calls"]:
                                read_set["calls"].append(sv_call.__name__)


blur_amount = 10
#blur_amount = 350

##
# id_b = ground truth
# id_a = calls
def analyze_by_score(call_table_analyzer, call_table, id_a, id_b):
    num_calls_a = call_table.num_calls(id_a, 0) # num calls made
    if num_calls_a == 0:
        return [], [], [], 0, 0
    min_score = call_table.min_score(id_a)
    max_score = call_table.max_score(id_a)
    inc = (max_score - min_score) / 50
    if max_score <= min_score:
        inc = 1
    p = max_score + 1
    #xs = []
    xs_2 = []
    #ys = []
    ys_2 = []
    ps = []
    num_calls_b = call_table.num_calls(id_b, 0) # num actual calls
    if num_calls_b == 0 or min_score == float('inf') or max_score == float('inf'):
        #print(num_calls_a, min_score, max_score)
        return [], [], [], 0, 0
    num_almost_overlaps_b_to_a = 0
    while p > min_score:
        #print(min_score, p, max_score)
        ps.append(max(p - inc, 0))
        # how many of the sv's are detected?
        #num_overlaps_b_to_a = call_table.num_overlaps(id_b, id_a, p, 0)
        num_almost_overlaps_b_to_a += call_table_analyzer.num_overlaps(id_a, id_b, p - inc, p, blur_amount)

        #xs.append(num_overlaps_b_to_a/num_calls_b)
        xs_2.append(num_almost_overlaps_b_to_a/num_calls_b)

        num_calls_a = call_table.num_calls(id_a, p - inc) # num calls made

        assert num_almost_overlaps_b_to_a <= num_calls_b
        assert num_almost_overlaps_b_to_a <= num_calls_a

        if num_calls_a == 0:
            #ys.append(0)
            ys_2.append(0)
        else:
            #ys.append(num_overlaps_b_to_a/num_calls_a)
            # assert(num_almost_overlaps_b_to_a <= num_calls_a)
            ys_2.append(num_almost_overlaps_b_to_a/num_calls_a)

        p -= inc

    # we have computed the lists in reverse order -> fix that
    xs_2.reverse()
    ys_2.reverse()
    ps.reverse()

    #num_invalid_calls = call_table.num_invalid_calls(id_a, 0, 0)
    num_invalid_calls_fuzzy = call_table_analyzer.num_invalid_calls(id_a, 0, max_score + 1, blur_amount)
    avg_blur = round(call_table_analyzer.blur_on_overlaps(id_a, id_b, 0, max_score + 1, blur_amount), 1)
    assert avg_blur <= blur_amount

    return xs_2, ys_2, ps, num_invalid_calls_fuzzy, avg_blur
    #return xs, ys, xs_2, ys_2, ps, num_invalid_calls, num_invalid_calls_fuzzy, avg_blur

##
# prints a bar plot that shows the ratio of false and true positives with respect the the diagonal shape filter of
# sv calls. 
# Diagonal shape filter:
#  - get standard deviation of sv jumps from mean delta value
#  - get standard deviation of sv jumps from mean r+q value
#  - divide values and filer calls based on that ratio
def compute_diagonal_threshold_picture(db_conn, id_a, id_b):
    parameter_set_manager = ParameterSetManager()
    true_positives = SvCallsFromDb(parameter_set_manager, db_conn, id_a, id_b, True, blur_amount)
    false_positives = SvCallsFromDb(parameter_set_manager, db_conn, id_a, id_b, False, blur_amount)

    def get_std(l):
        l.sort()
        c = len(l) // 2
        if len(l) % 2 == 1:
            mean = l[c]
        else:
            mean = (l[c-1] + l[c]) // 2
        sq_diff = 0
        for x in l:
            sq_diff += int((mean - x) ** 2)
        return sq_diff // len(l)

    def get_diagonal_value(call):
        diagonalA = []
        diagonalB = []
        for idx in range(len(call.supporing_jump_ids)):
            jump = call.get_jump(idx)
            x = jump.from_pos
            y = jump.to_pos
            diagonalA.append(y - x)
            diagonalB.append(y + x)
        return get_std(diagonalA) // max(get_std(diagonalB), 1)


    values_true = []
    while true_positives.hasNext():
        call = true_positives.next()
        values_true.append(get_diagonal_value(call))
    values_false = []
    while false_positives.hasNext():
        call = false_positives.next()
        #values_false.append(get_diagonal_value(call))
        values_false.append(0)

    values_true.sort()
    values_false.sort()

    num_buckets = 100
    min_ = min(values_true[0], values_false[0])
    max_ = max(values_true[-1], values_false[-1]) + 1
    div = max_ - min_
    def to_buckets(l):
        buckets = [0]*num_buckets
        for val in l:
            idx = num_buckets*(val - min_) // div
            #print(val, min_, max_, div, idx)
            assert idx < num_buckets
            buckets[idx] += 1
        return buckets
    
    buckets_true = to_buckets(values_true)
    buckets_false = to_buckets(values_false)

    plot = figure()
    plot.vbar(x=[x*div/num_buckets+min_ for x in range(num_buckets)], top=buckets_false, bottom=0,
              width=div*3/(4*num_buckets), legend="False positives", color="red", fill_alpha=0.5)
    plot.vbar(x=[x*div/num_buckets+min_ for x in range(num_buckets)], top=buckets_true, bottom=0,
              width=div*3/(4*num_buckets), legend="True positives", color="green", fill_alpha=0.5)
    plot.xaxis[0].formatter = PrintfTickFormatter(format="%d")
    show(plot)


def compare_caller(call_table_analyzer, call_table, id_a, id_b, min_score):
    num_calls_a = call_table.num_calls(id_a, min_score) # num calls made
    num_calls_b = call_table.num_calls(id_b, min_score) # num actual calls
    if num_calls_a == 0:
        print("no calls")
        return (0, 0, 0, 0, 0, 0, 0)
    max_score = call_table.max_score(id_a)
    call_area_a = call_table.call_area(id_a, min_score)
    if num_calls_a > 0:
        rel_call_area_a = int(math.sqrt(call_area_a/num_calls_a)) # get the edge length
    else:
        rel_call_area_a = "n/a"
    call_area_b = call_table.call_area(id_b, min_score)
    rel_call_area_b = call_area_b/num_calls_b
    # @note start with the largest blur_amount so that we only need to initialize the cache once
    num_almost_overlaps_a_to_b = call_table_analyzer.num_overlaps(id_a, id_b, min_score, max_score+1, blur_amount) # true positives
    # true positives
    num_overlaps_a_to_b = call_table_analyzer.num_overlaps(id_a, id_b, min_score, max_score+1, 0)
    num_almost_overlaps_a_to_b -= num_overlaps_a_to_b
    # how many of the sv's are detected?
    num_almost_overlaps_b_to_a = call_table_analyzer.num_overlaps(id_b, id_a, min_score, max_score+1, blur_amount)
    # how many of the sv's are detected?
    num_overlaps_b_to_a = call_table_analyzer.num_overlaps(id_b, id_a, min_score, max_score+1, 0)
    num_errors = num_calls_b - num_overlaps_b_to_a # how many of the sv's are NOT detected?
    num_way_off = num_calls_b - num_almost_overlaps_b_to_a # how many of the sv's are NOT detected?
    num_invalid_calls = call_table_analyzer.num_invalid_calls(id_a, min_score, max_score+1, 0)
    return (num_calls_a, num_overlaps_a_to_b, num_almost_overlaps_a_to_b, num_errors, num_way_off, 
            rel_call_area_a, num_invalid_calls)


def compare_all_callers_against(dataset_name, json_info_file, out_file_name=None, outfile_2_name=None):
    db_conn = DbConn(dataset_name)
    pool = PoolContainer(ParameterSetManager().get_num_threads() + 1, dataset_name)
    call_table_analyzer = SvCallTableAnalyzer(pool)
    run_table = SvCallerRunTable(db_conn)
    call_table = SvCallTable(db_conn)
    out = []
    out_2 = {}
    for dataset in json_info_file["datasets"]:
        id_b = dataset["ground_truth"]
        name_b = dataset["name"]
        date_b = run_table.getDate(id_b)
        print("ground truth is:", name_b, "-", date_b, "[ id:", id_b, "] - with",
              call_table.num_calls(id_b, 0), "calls")
        #print("sensitivity = true positive rate = recall")
        #print("missing rate = how many calls are missing")
        for id_a in run_table.newest_unique_runs(2, "ground_truth=" + str(id_b)):
            name_a = run_table.getName(id_a).split("-")
            dataset_name = name_a[0]
            dataset_size = name_a[1]
            seq = name_a[2]
            cov = name_a[3]
            aligner = name_a[4]
            caller = name_a[5]
            date_a = run_table.getDate(id_a)
            print("analyzing", id_a, name_a, date_a)
            out.append([dataset_name, dataset_size, seq, cov, caller, aligner, str(id_a),
                        *(str(x) for x in compare_caller(call_table_analyzer, call_table, id_a, id_b, 0))])
            if str(name_a[:4]) not in out_2:
                out_2[str(name_a[:4])] = {}
            out_2[str(name_a[:4])][str((aligner, caller))] = analyze_by_score(call_table_analyzer, call_table, id_a, id_b)
            #compute_diagonal_threshold_picture(db_conn, id_a, id_b)

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


def analyze_sample_dataset(dataset_name, run_callers=True, recompute_jumps=False, out_file_name=None, run_ma=True,
                           run_others=True, recompute_calls=False):
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
    with open(sv_data_dir + dataset_name + "/info.json", "r") as json_file:
        json_info_file = json.loads(json_file.read(), object_hook=_decode)

    # create the calls
    if run_callers:
        pack = Pack()
        pack.load(json_info_file["reference_path"] + "/ma/genome")
        fm_index = FMIndex()
        fm_index.load(json_info_file["reference_path"] + "/ma/genome")

        # create alignment files if they do not exist
        create_alignments_if_necessary(dataset_name, json_info_file, pack, fm_index, recompute_jumps, run_ma,
                                       run_others)
        # save the info.json file
        print(json_info_file)
        with open(sv_data_dir + dataset_name + "/info.json", "w") as json_out:
            json.dump(json_info_file, json_out)


        run_callers_if_necessary(dataset_name, json_info_file, pack, fm_index, run_others, recompute_calls)

        # save the info.json file
        print(json_info_file)
        with open(sv_data_dir + dataset_name + "/info.json", "w") as json_out:
            json.dump(json_info_file, json_out)

    compare_all_callers_against(dataset_name, json_info_file, sv_data_dir + dataset_name + "/bar_diagrams.tsv",
                                sv_data_dir + dataset_name + "/by_score.json")


if __name__ == "__main__":
    analyze_sample_dataset("comprehensive")