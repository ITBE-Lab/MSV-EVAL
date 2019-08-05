from bokeh.layouts import column
from bokeh.plotting import figure, show, reset_output, ColumnDataSource
from bokeh.models import Arrow, VeeHead
from bokeh.palettes import Category20, Category10
from create_json import create_json_from_db
import json
from MA import *
from sweep_sv_jumps import sv_jumps_to_dict

def render_from_dict(json_dict, start=None, end=None, on_y_aswell=True):
    # zip that just loops through the all shorter items, including scalars
    def zip_longest_scalar(*args):
        for i in range(max(len(x) if isinstance(x, list) else 1 for x in args)):
            yield (x[i%len(x)] if isinstance(x, list) else x for x in args)

    x_offset = json_dict["x_offset"]

    plots = []
    for panel in json_dict["panels"]:
        plot = figure(
            width=1700,
            height=panel["h"],
            tooltips="@i",
            #tools=[
            #    "pan", "xpan", "wheel_zoom", "xwheel_zoom", "box_zoom", "save",
            #    "reset", "hover"
            #],
            #active_drag="xpan",
            #active_scroll="xwheel_zoom"
            tools=[
                "pan", "wheel_zoom", "box_zoom", "save",
                "reset", "hover"
            ],
            active_drag="pan",
            active_scroll="wheel_zoom")
        if not len(plots) == 0:
            plot.x_range = plots[0].x_range
        elif not start is None and not end is None:
            plot.x_range.start = start
            plot.x_range.end = end
            if on_y_aswell:
                plot.y_range.start = start
                plot.y_range.end = end
        for item in panel["items"]:
            if item["type"] == "box":
                cds = {
                    'l': [],
                    'r': [],
                    't': [],
                    'b': [],
                }
                for x, y, w, h in item["data"]:
                    x += x_offset
                    cds["l"].append(x)
                    cds["r"].append(x + w)
                    cds["b"].append(y)
                    cds["t"].append(y + h)
                plot.quad(
                    left='l',
                    right='r',
                    top='t',
                    bottom='b',
                    line_width=0,
                    color=item['color'],
                    source=ColumnDataSource(cds))
            if item["type"] == "patch":
                cds = {
                    'xs': [],
                    'ys': []
                }
                for xs, ys in item["data"]:
                    cds["xs"].append([x + x_offset for x in xs])
                    cds["ys"].append(ys)
                plot.patches(
                    xs='xs',
                    ys='ys',
                    line_width=1,
                    color=item['color'],
                    source=ColumnDataSource(cds))
            elif item["type"] == "box-alpha":
                cds = {
                    'l': [],
                    'r': [],
                    't': [],
                    'b': [],
                    'a': [],
                    'i': []
                }
                for x, y, w, h, a, i in item["data"]:
                    x += x_offset
                    cds["l"].append(x)
                    cds["r"].append(x + w)
                    cds["b"].append(y)
                    cds["t"].append(y + h)
                    cds["a"].append(a)
                    cds["i"].append(i)
                plot.quad(
                    left='l',
                    right='r',
                    top='t',
                    bottom='b',
                    fill_alpha="a",
                    line_width=item['line_width'],
                    fill_color=item['color'],
                    line_color=item['line_color'],
                    source=ColumnDataSource(cds))
            elif item["type"] == "line":
                cds = {
                    "x": [],
                    "y": [],
                    'i': []
                }
                for i, (x, y, w, h) in enumerate(item["data"]):
                    x += x_offset
                    cds["x"].append([x, x+w])
                    cds["y"].append([y, y+h])
                    cds["i"].append(str(i))
                plot.multi_line(
                    xs='x',
                    ys='y',
                    line_width=3,
                    color=item['color'],
                    source=ColumnDataSource(cds))
            elif item["type"] == "plus":
                cds = {
                    "x": [],
                    "y": [],
                    "i": []
                }
                for x, y, i in item["data"]:
                    x += x_offset
                    cds["x"].append(x)
                    cds["y"].append(y)
                    cds["i"].append(i)
                plot.cross(
                    x='x',
                    y='y',
                    size=20,
                    line_width=3,
                    color=item['color'],
                    source=ColumnDataSource(cds))
            elif item["type"] == "arrow":
                for x, y, w, h in item["data"]:
                    x += x_offset
                    plot.add_layout(Arrow(end=VeeHead(size=10, fill_color=None, line_color=item["color"]),
                                                    x_start=x, y_start=y, x_end=x + w, y_end=y + h,
                                                    line_dash=[1, 2], line_width=3, line_color=item["color"]))
        plot.xaxis.visible = False
        plots.append(plot)

    plots[-1].xaxis.visible = True
    reset_output()
    show(column(plots))

def render_seeds(sv_db, seq_id, fm_index, out_dict, start, size):
    parameter_set_manager = ParameterSetManager()
    parameter_set_manager.set_selected("SV-Illumina")
    #parameter_set_manager.set_selected("SV-PacBio")
    #parameter_set_manager.set_selected("SV-ONT")

    lock_module = Lock(parameter_set_manager)
    seeding_module = BinarySeeding(parameter_set_manager)
    seeds_filter = libMA.FilterSeedsByArea(parameter_set_manager, start, size)
    collector = libMA.SeedsCollector(parameter_set_manager)

    fm_pledge = Pledge()
    fm_pledge.set(fm_index)

    res = VectorPledge()
    # graph for single reads
    for idx in range(parameter_set_manager.get_num_threads()):
        nuc_seq_getter = AllNucSeqFromSql(parameter_set_manager, sv_db, seq_id, idx,
                                          parameter_set_manager.get_num_threads())
        queries_pledge = promise_me(nuc_seq_getter)
        query_pledge = promise_me(lock_module, queries_pledge)
        segments_pledge = promise_me(seeding_module, fm_pledge, query_pledge)
        seeds_pledge = promise_me(seeds_filter, segments_pledge, fm_pledge, query_pledge)
        empty = promise_me(collector, seeds_pledge)
        unlock_pledge = promise_me(UnLock(parameter_set_manager, query_pledge), empty)
        res.append(unlock_pledge)

    # drain all sources
    res.simultaneous_get( parameter_set_manager.get_num_threads() )
    
    seeds = []
    curr_seeds = [s for s in collector.cpp_module.collection]
    curr_seeds.sort(key=lambda x: x.start)
    seeds_total = len(curr_seeds)
    print("rendering", seeds_total, "seeds")
    for idx, s in enumerate(curr_seeds):
        seeds.append((s, s.delta, idx))
    seeds.sort(key=lambda x: x[0].start_ref)
    ends_list = []

    render_list = [[] for _ in range(100)]

    for seed, nuc_seq_id, seed_idx in seeds:
        min_start_ref = seed.start_ref
        if not seed.on_forward_strand:
            min_start_ref = seed.start_ref - seed.size
        idx = len(ends_list)
        for index, end_pos in enumerate(ends_list):
            if end_pos + 3 < min_start_ref:
                idx = index
                break
        if idx == len(ends_list):
            ends_list.append(0)

        if seed.on_forward_strand:
            ends_list[idx] = seed.start_ref + seed.size
        else:
            ends_list[idx] = seed.start_ref
        
        xs = [min_start_ref, idx, # x, y
              seed.size, 0.5, # w, h
              1, # alpha
              "read_id: " + str(nuc_seq_id) + " q_pos: " + str(seed.start) + 
              " index on query: " + str(seed_idx) + 
              str(" forw" if seed.on_forward_strand else " rev")]
        render_list[nuc_seq_id%100].append(xs)


    panel = {
                "items": [
                    {
                        "type": "box-alpha",
                        "color": Category10[10][idx%10],
                        "line_color": Category10[10][int(idx/10)%10],
                        "line_width": 3,
                        "group": "seeds",
                        "data": render_list[idx]
                    } for idx in range(100)
                ],
                "h": 300
            }

    out_dict["panels"].append(panel)

    return out_dict


def render_from_json(json_file_name, start, end, on_y_aswell=False):
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
    with open(json_file_name, "r") as json_file:
        json_dict = json.loads(json_file.read(), object_hook=_decode)
    render_from_dict(json_dict, start, end, on_y_aswell)

if __name__ == "__main__":

    pos = 0# int(2.7835 * 10**9)
    size = int(10**10)# int(0.002 * 10**9)
    print(pos, size)

    #sv_db = SV_DB("/MAdata/databases/sv_simulated", "open")
    sv_db = SV_DB("/MAdata/sv_datasets/minimal-2/svs.db", "open")
    #out_dict = create_json_from_db(sv_db, "/MAdata/genome/human/GRCh38.p12/ma/genome")
    #out_dict = sv_jumps_to_dict(sv_db, [1, 12], pos, pos, size, size)
    out_dict = sv_jumps_to_dict(sv_db, [1, 2], only_supporting_jumps=True)

    #out_dict = sv_jumps_to_dict(sv_db, [1, 2], pos, pos, size, size, min_score=-1)

    #fm_index = FMIndex()
    #fm_index.load("/MAdata/genome/random_w_mobile_ele/ma/genome")
    #fm_index.load("/MAdata/genome/random_10_pow_6/ma/genome")
    #fm_index.load("/MAdata/genome/human/GRCh38.p12/ma/genome")
    #out_dict = render_seeds(sv_db, 1, fm_index, out_dict, pos, size)

    render_from_dict(out_dict)