from bokeh.layouts import column
from bokeh.plotting import figure, show, reset_output, ColumnDataSource
from bokeh.models import Arrow, VeeHead
import json

def render_from_json(json_file_name, start, end):
    # zip that just loops through the all shorter items, including scalars
    def zip_longest_scalar(*args):
        for i in range(max(len(x) if isinstance(x, list) else 1 for x in args)):
            yield (x[i%len(x)] if isinstance(x, list) else x for x in args)
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

    # first add x_offset to all x values (to keep python code simple)
    x_offset = json_dict["x_offset"]
    for panel in json_dict["panels"]:
        for item in panel["items"]:
            item["x"] = [x + x_offset for x in item["x"]]

    plots = []
    for panel in json_dict["panels"]:
        plot = figure(
            width=1700,
            height=panel["h"],
            x_range=[start, end] if len(plots) == 0 else plots[0].x_range,
            tools=[
                "pan", "xpan", "wheel_zoom", "xwheel_zoom", "box_zoom", "save",
                "reset"
            ],
            active_drag="xpan",
            active_scroll="xwheel_zoom")
        for item in panel["items"]:
            if item["type"] == "box":
                cds = {
                    'l': [],
                    'r': [],
                    't': [],
                    'b': [],
                }
                for x, y, w, h in zip_longest_scalar(item["x"], item["y"], item["w"], item["h"]):
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
            elif item["type"] == "line":
                cds = {
                    "x": [],
                    "y": []
                }
                for x, y, w, h in zip_longest_scalar(item["x"], item["y"], item["w"], item["h"]):
                    cds["x"].append([x, x+w])
                    cds["y"].append([y, y+h])
                plot.multi_line(
                    xs='x',
                    ys='y',
                    line_width=3,
                    color=item['color'],
                    source=ColumnDataSource(cds))
            elif item["type"] == "arrow":
                for x, y, w, h in zip_longest_scalar(item["x"], item["y"], item["w"], item["h"]):
                    plot.add_layout(Arrow(end=VeeHead(size=10, fill_color=None, line_color=item["color"]),
                                                    x_start=x, y_start=y, x_end=x + w, y_end=y + h,
                                                    line_dash=[1, 2], line_width=3, line_color=item["color"]))
        plot.xaxis.visible = False
        plot.yaxis.visible = False
        plots.append(plot)

    plots[-1].xaxis.visible = True
    reset_output()
    show(column(plots))