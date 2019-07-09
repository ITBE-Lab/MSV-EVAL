from bokeh.layouts import column
from bokeh.plotting import figure, show, reset_output, ColumnDataSource
from bokeh.models import Arrow, VeeHead
from create_json import create_json_from_db
import json
from MA import SV_DB, ParameterSetManager
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
    #sv_db = SV_DB("/MAdata/databases/sv_simulated", "open")
    sv_db = SV_DB("/MAdata/sv_datasets/minimal/svs.db", "open")
    #out_dict = create_json_from_db(sv_db, "/MAdata/genome/human/GRCh38.p12/ma/genome")
    out_dict = sv_jumps_to_dict(sv_db, [1, 32], 600000, 600000, 50000, 50000)
    
    #out_dict = sv_jumps_to_dict(sv_db, [3, 122], only_supporting_jumps=True)
    render_from_dict(out_dict)