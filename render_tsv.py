from bokeh.layouts import layout
from bokeh.plotting import figure, show, reset_output, ColumnDataSource
from bokeh.models import Arrow, VeeHead, FactorRange, LabelSet
from bokeh.palettes import Category20, Category10
from bokeh.transform import dodge
from bokeh.core.properties import value
from create_json import create_json_from_db
import json
from MA import SV_DB, ParameterSetManager
from sweep_sv_jumps import sv_jumps_to_dict
from bokeh.models import Legend

def split_by_cat(s, e, l):
    a = 0
    b = 0
    while a < len(l):
        while b < len(l) and l[b][s:e+1] == l[a][s:e+1]:
            b += 1
        yield l[a][s:e+1], l[a:b]
        a = b

def to_int(x):
    if x == "n/a":
        return 0
    return int(x)

def render_from_list(tsv_list, json_dict, plot_category=(0,0), plot_sub_category=(0,2), category=(3,5), 
                     stacks=[(8, 9, 11, "green", "orange", "gray")],
                     bars=[(7, 9, "red"), (12, None, "blue"), (13, None, "purple")]):

    plotss = []
    for _, sub_lists in split_by_cat(*plot_category, tsv_list[1:]):
        plots = []
        for name_, sub_list in split_by_cat(*plot_sub_category, sub_lists):
            name = str(name_)
            x = []
            for row in sub_list:
                tup = (*row[category[0] : category[1]], row[category[1]] + " [id: " + str(row[6]) + "]")
                x.append(tup)
            range_num = to_int(sub_list[0][8])+to_int(sub_list[0][10])
            plot = figure(x_range=FactorRange(*x),
                          y_range=(-range_num*0.1, range_num * 1.1),
                          title=name,
                          tools="save",
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

                # helper lines
                #plot.line(x=[x[0], x[-1]], 
                #          y=[to_int(sub_list[0][bottom2_idx]), to_int(sub_list[0][bottom2_idx])], 
                #          color="orange", line_dash=[3,1])
                #plot.line(x=[x[0], x[-1]], 
                #          y=[to_int(sub_list[0][bottom_idx]), to_int(sub_list[0][bottom_idx])], 
                #          color="green", line_dash=[3,1])

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

    for name_1, sub_lists in split_by_cat(0, 1, tsv_list[1:]):
        plotss = []
        for name_2, sub_lists_2 in split_by_cat(2, 2, sub_lists):
            plotss.append([])
            for name_3, sub_list in split_by_cat(3, 3, sub_lists_2):
                legend = []
                name = str([*name_1, *name_2, *name_3])
                plot_3 = figure(title=name + " - 10nt blur", tooltips="@i", active_drag=None, width=500, height=300)
                for idx, row in enumerate(sub_list):
                    x_every_2 = 1
                    x_every = 4
                    aligner_name = str((row[5], row[4]))
                    if aligner_name in json_dict[name]:
                        x_2, y_2, p, num_invalid_calls_fuzzy, avg_blur = json_dict[name][aligner_name]

                        if len(x_2) > x_every_2:
                            line = plot_3.line(x="x", y="y", color=Category10[10][idx%10],
                                        source=ColumnDataSource(data=dict(x=x_2[::x_every_2], y=y_2[::x_every_2], 
                                                                        i=p[::x_every_2])),
                                        line_width=3, alpha=0.5)
                        else:
                            line = plot_3.x(x="x", y="y", color=Category10[10][idx%10],
                                        source=ColumnDataSource(data=dict(x=x_2[::x_every_2], y=y_2[::x_every_2], 
                                                                        i=p[::x_every_2])),
                                        line_width=3, alpha=0.5, size=8)
                        legend.append( (aligner_name + " #inv:" + str(num_invalid_calls_fuzzy) + " #avgdist: " + \
                                        str(avg_blur), [line]) )
                if len(plotss[0]) > 0:
                    plot_3.x_range = plotss[0][0].x_range
                    plot_3.y_range = plotss[0][0].y_range
                    
                plot_3.xaxis.axis_label = "recall"
                plot_3.yaxis.axis_label = "precision"
                legend_obj = Legend(items=legend, location="center")
                plot_3.add_layout(legend_obj, "right")
                plotss[-1].append(plot_3)

        #for name_2, sub_lists_2 in split_by_cat(2, 2, sub_lists):
        #    plotss.append([])
        #    for name_3, sub_list in split_by_cat(3, 3, sub_lists_2):
        #        legend = []
        #        name = str([*name_1, *name_2, *name_3])
        #        plot_2 = figure(title=name, tooltips="@i", active_drag=None, width=500, height=300)
        #        for idx, row in enumerate(sub_list):
        #            x_every_2 = 1
        #            x_every = 4
        #            aligner_name = str((row[5], row[4]))
        #            if aligner_name in json_dict[name]:
        #                x, y, x_2, y_2, p, num_invalid_calls, _ = json_dict[name][aligner_name]
        #
        #                if len(x) > x_every_2:
        #                    line = plot_2.line(x="x", y="y", color=Category10[10][idx%10],
        #                                source=ColumnDataSource(data=dict(x=x[::x_every_2], y=y[::x_every_2], 
        #                                                        i=p[::x_every_2])),
        #                                line_width=3, alpha=0.5)
        #                else:
        #                    line = plot_2.x(x="x", y="y", color=Category10[10][idx%10],
        #                                source=ColumnDataSource(data=dict(x=x[::x_every_2], y=y[::x_every_2], 
        #                                                        i=p[::x_every_2])),
        #                                line_width=3, alpha=0.5, size=8)
        #                legend.append((aligner_name + " #inv:" + str(num_invalid_calls), [line]))
        #        if len(plotss[0]) > 0:
        #            plot_2.x_range = plotss[0][0].x_range
        #            plot_2.y_range = plotss[0][0].y_range
        #        plot_2.xaxis.axis_label = "recall"
        #        plot_2.yaxis.axis_label = "precision"
        #        legend_obj = Legend(items=legend, location="center")
        #        plot_2.add_layout(legend_obj, "right")
        #        plotss[-1].append(plot_2)

        reset_output()
        show(layout(plotss))

def print_ground_truth(dataset_name):
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
    sv_db = SV_DB("/MAdata/sv_datasets/" + dataset_name + "/svs.db", "open")
    for dataset in json_info_file["datasets"]:
        id_b = dataset["ground_truth"]
        name_b = dataset["name"]
        date_b = sv_db.get_run_date(id_b)
        print("ground truth is:", name_b, "-", date_b, "[ id:", id_b, "] - with", sv_db.get_num_calls(id_b, 0), "calls")

def render_from_tsv(dataset_name):
    print_ground_truth(dataset_name)

    tsv_list = []
    with open("/MAdata/sv_datasets/" + dataset_name + "/bar_diagrams.tsv", "r") as tsv_file:
        for line in tsv_file:
            tsv_list.append(line.split("\t"))
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
    #actually open and load the file
    json_dict = None # noop
    with open("/MAdata/sv_datasets/" + dataset_name + "/by_score.json", "r") as json_file:
        json_dict = json.loads(json_file.read(), object_hook=_decode)
    fist_line = tsv_list[0]
    tsv_list = tsv_list[1:]
    tsv_list.sort()
    tsv_list.insert(0, fist_line)
    #print(tsv_list)
    render_from_list(tsv_list, json_dict)

if __name__ == "__main__":
    render_from_tsv("minimal")