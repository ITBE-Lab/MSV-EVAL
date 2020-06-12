from MSV import *
from bokeh.models import ColumnDataSource, LabelSet
from bokeh.plotting import figure, show
from bokeh.models import FactorRange
from bokeh.layouts import column, gridplot

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

def seed_plot(seeds, x_from, x_to, section_size, dataset_name, sections):
    seed_plot = figure(title="Seed Representation - " + dataset_name, x_range=(x_from, x_to), y_range=(x_from, x_to))
    forw_x = []
    forw_y = []
    rev_x = []
    rev_y = []
    for seed in seeds:
        if seed.on_forward_strand:
            forw_x.append(seed.start_ref)
            forw_x.append(seed.start_ref + seed.size)
            forw_x.append(float("NaN"))
            forw_y.append(seed.start)
            forw_y.append(seed.start + seed.size)
            forw_y.append(float("NaN"))
        else:
            rev_x.append(seed.start_ref)
            rev_x.append(seed.start_ref - seed.size)
            rev_x.append(float("NaN"))
            rev_y.append(seed.start)
            rev_y.append(seed.start + seed.size)
            rev_y.append(float("NaN"))
    seed_plot.line(x=forw_x, y=forw_y, line_width=1)
    seed_plot.line(x=rev_x, y=rev_y, line_color="orange", line_width=1)
    seed_plot.xaxis.ticker = sections
    seed_plot.yaxis.ticker = sections
    seed_plot.xgrid.ticker = seed_plot.xaxis.ticker
    seed_plot.ygrid.ticker = seed_plot.yaxis.ticker
    return seed_plot

def vcf_pic(x_range, from_to_calls_lists, dataset_name, sections):
    ys = []
    lefts = []
    center = []
    rights = []
    names = []
    ids = []
    for name, from_to_calls_list in from_to_calls_lists:
        for from_pos, to_pos, call_id, call_name in from_to_calls_list:
            ys.append((name, str(call_id)))
            lefts.append(from_pos)
            center.append((from_pos + to_pos) / 2)
            rights.append(to_pos)
            names.append(call_name)

    plot = figure(title='VCF Output - ' + dataset_name, x_range=x_range, y_range=FactorRange(*sorted(ys,reverse=True)))
    plot.hbar(y=ys, left=lefts, right=rights, height=.8)
    plot.add_layout(LabelSet(x='c', y='y', text='n', level='annotation', text_align='center', text_baseline='middle',
                             source=ColumnDataSource(data={'c':center, 'y':ys, 'n':names}), render_mode='canvas'))
    plot.xaxis.ticker = sections
    plot.xgrid.ticker = plot.xaxis.ticker
    return plot

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
    detail_plots1 = []
    detail_plots2 = []
    for caller_ids, read_types, gt_id, seeds, from_to_calls_lists, x_from, x_to, section_size in sets:
        sections = [*range(x_from, x_to+1, section_size)]
        dataset_name = run_table.getName(gt_id)
        detail_plots1.append( seed_plot(seeds, x_from, x_to, section_size, dataset_name, sections) )
        detail_plots2.append(vcf_pic(detail_plots1[-1].x_range, from_to_calls_lists, dataset_name, sections))
        repl_name_w_color = {
            "match": "green",
            "missing": "red",
            "one_axis_match": "orange",
            "additional": "magenta"
        }
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
    show(column(plot, gridplot([detail_plots1, detail_plots2], plot_width=400, plot_height=400)))




