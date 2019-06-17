from bokeh.layouts import layout
from bokeh.plotting import figure, show, reset_output, ColumnDataSource
from bokeh.models import Arrow, VeeHead, FactorRange, LabelSet
from bokeh.transform import dodge
from bokeh.core.properties import value
from create_json import create_json_from_db
import json
from MA import SV_DB, ParameterSetManager
from sweep_sv_jumps import sv_jumps_to_dict


def split_by_cat(s, e, l):
    a = 0
    b = 0
    while a < len(l):
        while b < len(l) and l[b][s:e+1] == l[a][s:e+1]:
            b += 1
        yield str(l[a][s:e+1]), l[a:b]
        a = b

def to_int(x):
    if x == "n/a":
        return 0
    return int(x)

def render_from_list(tsv_list, plot_category=(0,0), plot_sub_category=(0,2), category=(3,5), 
                     stacks=[(8, 9, 11, "green", "orange", "gray")],
                     bars=[(7, 9, "red"), (12, None, "blue")]):

    plotss = []
    for _, sub_lists in split_by_cat(*plot_category, tsv_list[1:]):
        plots = []
        for name, sub_list in split_by_cat(*plot_sub_category, sub_lists[1:]):
            x = []
            for row in sub_list:
                tup = (*row[category[0] : category[1]], row[category[1]] + " [id: " + str(row[6]) + "]")
                x.append(tup)
            range_num = to_int(sub_list[0][8])+to_int(sub_list[0][10])
            plot = figure(x_range=FactorRange(*x),
                          y_range=(-range_num*0.1, range_num * 1.1),
                          title=name,
                          toolbar_location=None,
                          tools="",
                          width=800)
            n = len(stacks) + len(bars)
            for idx, (bottom_idx, bottom2_idx, top_idx, bottom_color, bottom2_color, top_color) in enumerate(stacks):
                data = { 'x':x,
                        "bottom":[to_int(row[bottom_idx]) for row in sub_list],
                        "bottom2":[to_int(row[bottom2_idx]) - to_int(row[bottom_idx]) for row in sub_list],
                        "top":[to_int(row[top_idx]) for row in sub_list],

                        "a":[to_int(row[bottom2_idx]) + to_int(row[top_idx])/2 if to_int(row[top_idx]) != 0 else -100 for row in sub_list],
                        "b":[to_int(row[bottom_idx])/2 if to_int(row[bottom_idx]) != 0 else -100 for row in sub_list],
                        "b2":[to_int(row[bottom_idx]) + (to_int(row[bottom2_idx]) - to_int(row[bottom_idx]))/2 if row[bottom2_idx] != row[bottom_idx] else -100 for row in sub_list]
                    }
                plot.vbar_stack(["bottom", "bottom2", "top"],
                                x=dodge('x', -0.4 + idx / n, range=plot.x_range),
                                width=0.8/n,
                                color=[bottom_color, bottom2_color, top_color],
                                source=ColumnDataSource(data=data))#,
                                #legend=[value(tsv_list[0][bottom_idx]), value(tsv_list[0][bottom2_idx]),
                                #        value(tsv_list[0][top_idx])])
                
                labels = LabelSet(x=dodge('x', -0.4 + idx / n, range=plot.x_range), y='a', text='top', level='glyph',
                                x_offset=0, y_offset=0, source=ColumnDataSource(data=data), render_mode='css',
                                text_align="center", text_color="white", angle=1.5708)
                plot.add_layout(labels)
                labels = LabelSet(x=dodge('x', -0.4 + idx / n, range=plot.x_range), y='b', text='bottom', level='glyph',
                                x_offset=0, y_offset=0, source=ColumnDataSource(data=data), render_mode='css',
                                text_align="center", text_color="white", angle=1.5708)
                plot.add_layout(labels)
                labels = LabelSet(x=dodge('x', -0.4 + idx / n, range=plot.x_range), y='b2', text='bottom2',
                                level='glyph', x_offset=0, y_offset=0, source=ColumnDataSource(data=data),
                                render_mode='css', text_align="center", text_color="white", angle=1.5708)
                plot.add_layout(labels)

            for idx, (data_idx, sub_idx, color) in enumerate(bars):
                name = tsv_list[0][data_idx]
                data = { 'x':x }
                if not sub_idx is None:
                    data[name] = [to_int(row[data_idx]) - to_int(row[sub_idx]) for row in sub_list]
                    data['a'] = [(to_int(row[data_idx]) - to_int(row[sub_idx]))/2 for row in sub_list]
                else:
                    data[name] = [to_int(row[data_idx]) for row in sub_list]
                    data['a'] = [to_int(row[data_idx])/2 for row in sub_list]

                plot.vbar(x=dodge('x', -0.4 + ( idx + len(stacks) ) / n, range=plot.x_range),
                        top=name,
                        width=0.8/n,
                        color=color,
                        source=ColumnDataSource(data=data))#,
                        #legend=value(name if sub_idx is None else "extra " + name))
                labels = LabelSet(x=dodge('x', -0.4 + ( idx + len(stacks) ) / n, range=plot.x_range), y='a',
                                  text=name, level='glyph', x_offset=0, y_offset=0, source=ColumnDataSource(data=data),
                                  render_mode='css', text_align="center", text_color="white", angle=1.5708)
                plot.add_layout(labels)


            plot.x_range.range_padding = 0.1
            #plot.y_range.range_padding = 0.1
            plot.xaxis.major_label_orientation = 1
            plot.xgrid.grid_line_color = None
            plot.legend.location = "top_center"
            plot.legend.orientation = "horizontal"
            plots.append(plot)
        plotss.append(plots)

    reset_output()
    show(layout(plotss))

def render_from_tsv(dataset_name):
    tsv_list = []
    with open("/MAdata/sv_datasets/" + dataset_name + "/bar_diagrams.tsv", "r") as tsv_file:
        for line in tsv_file:
            tsv_list.append(line.split("\t"))
    tsv_list.sort()
    render_from_list(tsv_list)

if __name__ == "__main__":
    render_from_tsv("comprehensive_random")