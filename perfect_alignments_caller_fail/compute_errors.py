from MSV import *
from bokeh.models import ColumnDataSource
from bokeh.plotting import figure, show
from bokeh.models import FactorRange

def call_to_points(db_con, caller_id):
    call_getter = SvCallsFromDb(ParameterSetManager(), db_con, caller_id)
    ret = []
    while call_getter.hasNext():
        call = call_getter.next()
        x = call.x.start + call.x.size/2
        y = call.y.start + call.y.size/2
        f = call.from_forward
        t = call.to_forward
        ret.append((x, y, f, t))
    return ret

def match(a, b, max_dist, num_axis=2):
    score = 0
    for idx in [0, 1]:
        if abs(a[idx] - b[idx]) < max_dist:
            score += 1
    return score >= num_axis and a[2:] == b[2:]

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


def render(db_con, sets, max_dist=10):
    x_range = []
    y_range = set()
    x = []
    y = []
    color = []
    label = []
    run_table = SvCallerRunTable(db_con)
    for caller_ids, read_types, gt_id in sets:
        repl_name_w_color = {
            "match": "green",
            "missing": "red",
            "one_axis_match": "orange",
            "additional": "magenta"
        }
        dataset_name = run_table.getName(gt_id)
        max_len = 0
        for caller_id, read_type in zip(caller_ids, read_types):
            caller_name = run_table.getName(caller_id)
            y_range.add((read_type, caller_name))
            for idx, error in enumerate(get_errors(db_con, caller_id, gt_id, max_dist)):
                if idx >= max_len:
                    max_len = idx + 1
                x.append((dataset_name, str(idx)))
                y.append((read_type, caller_name))
                color.append(repl_name_w_color[error])
                label.append(error)
        x_range.extend( [ (dataset_name, str(idx)) for idx in range(max_len)] )

    y_range.sort()

    plot = figure(title="Matching Breakpoints", x_range=FactorRange(*x_range), y_range=FactorRange(*y_range),
                  plot_width=800)
    plot.xaxis.axis_label = "Datasets & Break Points"
    plot.yaxis.axis_label = "Callers"
    r = plot.rect(x="x", y="y", width=0.8, height=0.8, line_color=None, legend_field='label', fill_color="color",
                  source=ColumnDataSource({"x":x, "y":y, "color":color, "label":label}))
    show(plot)




