from MSV import *
from bokeh.plotting import figure, show
from bokeh.models import FactorRange

def call_to_points(db_con, caller_id):
    call_getter = SvCallsFromDb(ParameterSetManager(), db_con, caller_id)
    ret = []
    while call_getter.hasNext():
        call = call_getter.next()
        x = call.x.start + call.x.size/2
        y = call.y.start + call.y.size/2
        s = call.switch_strand
        ret.append((x, y, s))
    return ret

def match(a, b, max_dist, num_axis=2):
    score = 0
    for idx in [0, 1]:
        if abs(a[idx] - b[idx]) < max_dist:
            score += 1
    return score >= num_axis

def find(l, c):
    for idx, ele in enumerate(l):
        if c(ele):
            return idx
    return None

def get_errors(db_con, caller_id, gt_id, max_dist=10):
    calls = call_to_points(db_con, caller_id)
    gts = call_to_points(db_con, gt_id)

    classification = ["missing"]*len(gts)

    # look for proper matches
    for i, gt in enumerate(gts):
        j = find(calls, lambda call: match(gt, call, max_dist))
        if not j is None:
            classification[i] = "match"
            del calls[j]

    # look for one axis matches
    for i, gt in enumerate(gts):
        if classification[i] == "missing":
            j = find(calls, lambda call: match(gt, call, max_dist, num_axis=1))
            if not j is None:
                classification[i] = "one_axis_match"
                del calls[j]

    # append left over calls
    for _ in calls:
        classification.append("additional")
    
    return classification


def render(db_con, caller_ids, gt_id, max_dist=10):
    classifications = []
    for caller_id in caller_ids:
        classifications.append(get_errors(db_con, caller_id, gt_id, max_dist))
    
    # make classifications same length
    max_len = max(len(x) for x in classifications)
    for idx in range(len(classifications)):
        while len(classifications[idx]) < max_len:
            classifications[idx].append("none")

    x_range = [str(idx) for idx in range(len(classifications[0]))]
    run_table = SvCallerRunTable(db_con)
    y_range = [run_table.getName(caller_id) for caller_id in caller_ids]

    repl_name_w_color = {
        "match": "green",
        "missing": "red",
        "one_axis_match": "orange",
        "additional": "magenta",
        "none": "white"
    }
    colors = [repl_name_w_color[item] for sublist in classifications for item in sublist]
    x = []
    y = []
    for y_ele in y_range:
        x.extend(x_range)
        y.extend([y_ele]*len(x_range))

    plot = figure(title=run_table.getName(gt_id), x_range=FactorRange(*x_range), y_range=FactorRange(*y_range),
                  plot_width=800)
    plot.xaxis.axis_label = "Break Points"
    plot.yaxis.axis_label = "Callers"
    plot.rect(x=x, y=y, width=0.8, height=0.8, line_color=None, fill_color=colors)
    show(plot)




