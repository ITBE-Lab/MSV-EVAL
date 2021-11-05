from MSV import *
from bokeh.models import ColumnDataSource, LabelSet
from bokeh.plotting import figure, show, save
from bokeh.models import FactorRange
from bokeh.layouts import column, gridplot
from bokeh.models.tools import HoverTool
from sv_util.settings import *

def call_to_points(x_range, entries, dataset_name, sections):
    plot = figure(title='Entries - ' + dataset_name, x_range=x_range, y_range=x_range)
    xs = []
    ys = []
    cs = []
    fs = []
    ts = []
    
    colors = {True:{True:"blue", False:"green"}, False:{True:"purple", False:"orange"}}
    for jump in entries:
        xs.append(jump.from_pos)
        ys.append(jump.to_pos)
        cs.append(colors[jump.from_forward][jump.to_forward])
        fs.append("fow" if jump.from_forward else "rev")
        ts.append("fow" if jump.to_forward else "rev")
    plot.x(x="xs", y="ys", color="cs", line_width=4, source=ColumnDataSource(data={
                    "xs":xs, "ys":ys, "cs":cs, "fs":fs, "ts":ts,
                }))
    plot.add_tools(HoverTool(tooltips=[("x, y", "@xs, @ys"), ("strandinfo", "from @fs to @ts")]))
    return plot

def seed_plot(seeds, dataset_name, sections):
    seed_plot = figure(title="Seed Representation - " + dataset_name,
                       x_range=(seeds[0].start_ref, seeds[-1].start_ref + seeds[-1].size))
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
    seed_plot.xgrid.ticker = seed_plot.xaxis.ticker
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

def render(sets, section_size):
    idx_to_abc = "abcdefghijklmnopqrstuvwxyz"
    x_range = []
    y_range = set()
    x = []
    y = []
    color = []
    label = []
    label_t = []
    detail_plots1 = []
    detail_plots2 = []
    detail_plots3 = []
    for seeds, from_to_calls_lists, dataset_name, msv_entries in sets:
        sections = [*range(seeds[0].start_ref, seeds[-1].start_ref + seeds[-1].size + 1, section_size)]
        detail_plots1.append(seed_plot(seeds, dataset_name, sections) )
        detail_plots2.append(vcf_pic(detail_plots1[-1].x_range, from_to_calls_lists, dataset_name, sections))
        detail_plots3.append(call_to_points(detail_plots1[-1].x_range, msv_entries, dataset_name, sections))
    if show_plots:
        show(gridplot([detail_plots1, detail_plots2, detail_plots3], plot_width=400, plot_height=400))
    if save_plots:
        save(gridplot([detail_plots1, detail_plots2, detail_plots3], plot_width=400, plot_height=400))


