from MA import *

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
        rel_call_area_a = call_area_a/num_calls_a
    else:
        rel_call_area_a = 0
    call_area_b = sv_db.get_call_area(id_b, min_score)
    rel_call_area_b = call_area_b/num_calls_b
    num_overlaps_a_to_b = sv_db.get_num_overlaps_between_calls(id_a, id_b, min_score, 0) # true positives
    num_almost_overlaps_a_to_b = sv_db.get_num_overlaps_between_calls(id_a, id_b, min_score, 100) # true positives
    num_almost_overlaps_a_to_b -= num_overlaps_a_to_b
    num_overlaps_b_to_a = sv_db.get_num_overlaps_between_calls(id_b, id_a, min_score, 0) # how many of the sv's are detected?
    num_errors = num_calls_b - num_overlaps_b_to_a # how many of the sv's are NOT detected?
    return (num_calls_a, num_overlaps_a_to_b, num_almost_overlaps_a_to_b, num_errors, rel_call_area_a)

def compare_callers(db_name, names_a, names_b=["simulated sv"], min_scores=[0]):
    sv_db = SV_DB(db_name, "open")
    #print("sensitivity = true positive rate = recall")
    #print("missing rate = how many calls are missing")
    print("test set", "ground truth", "min score", "#calls", "#found", "#almost", "#missed", "call-area", sep="\t")
    for name_a, name_b in zip(names_a, names_b):
        id_a = sv_db.get_run_id(name_a)
        id_b = sv_db.get_run_id(name_b)
        date_a = sv_db.get_run_date(id_a)
        date_b = sv_db.get_run_date(id_b)
        for min_score in min_scores:
            print(name_a + " - " + date_a, name_b + " - " + date_b, min_score,
                  *compare_caller(sv_db, id_a, id_b, min_score), sep="\t")

def compare_all_callers_against(db_name, name_b="simulated sv", min_scores=[0]):
    sv_db = SV_DB(db_name, "open")
    id_b = sv_db.get_run_id(name_b)
    date_b = sv_db.get_run_date(id_b)
    #print("sensitivity = true positive rate = recall")
    #print("missing rate = how many calls are missing")
    print("ground truth is: ", name_b, "-", date_b)
    out = [["test set", "time", "t", "#calls", "#found", "#almost", "#missed", "call-area"]]
    id_a = 0
    keep_looping = sv_db.get_num_runs()
    while keep_looping > 1:
        id_a += 1
        if not sv_db.run_exists(id_a):
            continue
        if id_a == id_b:
            continue
        name_a = sv_db.get_run_name(id_a)
        date_a = sv_db.get_run_date(id_a)
        for min_score in min_scores:
            out.append([name_a, date_a, str(min_score), 
                        *(str(x) for x in compare_caller(sv_db, id_a, id_b, min_score))])
        keep_looping -= 1
    print_columns(out)


#print("===============")
#compare_callers("/MAdata/databases/sv_simulated", ["MA-SV"])
#print("===============")
if __name__ == "__main__":
    compare_all_callers_against("/MAdata/databases/sv_simulated", min_scores=[10])