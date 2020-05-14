from MSV import *
from bokeh.models import ColumnDataSource, LabelSet
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
        ret.append((call.id, x, y, f, t))
    ret.sort()
    return ret

def match(a, b, max_dist, num_axis=2):
    score = 0
    for idx in [1, 2]:
        if abs(a[idx] - b[idx]) < max_dist:
            score += 1
    return score >= num_axis and a[3:] == b[3:]

def find(l, c):
    for idx, ele in enumerate(l):
        if c(ele):
            return idx
    return None

def get_errors(db_con, caller_id, gt_id, desc_table, max_dist=10):
    calls = call_to_points(db_con, caller_id)
    gts = call_to_points(db_con, gt_id)

    classification = [("missing", "")]*len(gts)

    # look for proper matches
    for i, gt in enumerate(gts):
        j = find(calls, lambda call: match(gt, call, max_dist))
        if not j is None:
            classification[i] = ("match", desc_table.get_desc(calls[j][0]) )
            del calls[j]

    # look for one axis matches
    for i, gt in enumerate(gts):
        if classification[i] == "missing":
            j = find(calls, lambda call: match(gt, call, max_dist, num_axis=1))
            if not j is None:
                classification[i] = ("one_axis_match", desc_table.get_desc(calls[j][0]) )
                del calls[j]

    # append left over calls
    for call in calls:
        classification.append( ("additional", desc_table.get_desc(call[0])) )

    return classification


def render(db_con, sets, max_dist=10):
    idx_to_abc = "abcdefghijklmnopqrstuvwxyz"
    x_range = []
    y_range = set()
    x = []
    y = []
    color = []
    label = []
    run_table = SvCallerRunTable(db_con)
    desc_table = CallDescTable(db_con)
    label_t = []
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
            if run_table.exists(caller_id):
                caller_name = run_table.getName(caller_id)
                y_range.add((read_type, caller_name))
                for idx, (error, desc) in enumerate(get_errors(db_con, caller_id, gt_id, desc_table, max_dist)):
                    if idx >= max_len:
                        max_len = idx + 1
                    x.append((dataset_name, idx_to_abc[idx]))
                    y.append((read_type, caller_name))
                    color.append(repl_name_w_color[error])
                    label.append(error)
                    label_t.append(desc)
        x_range.extend( [ (dataset_name, idx_to_abc[idx]) for idx in range(max_len)] )

    y_range = [*y_range]
    y_range.sort()

    plot = figure(title="Matching Breakpoints", x_range=FactorRange(*x_range), y_range=FactorRange(*y_range),
                  plot_width=800)
    plot.xaxis.axis_label = "Datasets & Break Points"
    plot.yaxis.axis_label = "Callers"
    plot.xaxis.group_label_orientation = math.pi/2
    plot.xaxis.subgroup_label_orientation = math.pi/2
    r = plot.rect(x="x", y="y", width=0.8, height=0.8, line_color=None, legend_field='label', fill_color="color",
                  source=ColumnDataSource({"x":x, "y":y, "color":color, "label":label}))
    labels = LabelSet(x='x', y='y', text='t', level='overlay', angle=math.pi/2, text_align="center",
              text_baseline="middle", x_offset=5, y_offset=5, source=ColumnDataSource({"x":x, "y":y, "t":label_t}),
              render_mode='canvas')
    plot.add_layout(labels)
    show(plot)




