# myapp.py

from random import random

from bokeh.layouts import column
from bokeh.models import Button
from bokeh.palettes import RdYlBu3
from bokeh.plotting import figure, curdoc
from bokeh.models.callbacks import CustomJS

args = curdoc().session_context.request.arguments
try:
    print(args.get("xs"))
    xs = float(args.get("xs")[0])
    ys = float(args.get("ys")[0])
    xe = float(args.get("xe")[0])
    ye = float(args.get("ye")[0])
except:
    xs = 0
    ys = 0
    xe = 100
    ye = 100


print(xs, xe, ys, ye)
# create a plot and style its properties
p = figure(x_range=(xs, xe), y_range=(ys, ye), toolbar_location=None)

# add a text renderer to our plot (no data yet)
r = p.text(x=[], y=[], text=[], text_color=[], text_font_size="20pt",
           text_baseline="middle", text_align="center")

i = 0

ds = r.data_source

new_data = dict()
new_data['x'] = ds.data['x'] + [10]
new_data['y'] = ds.data['y'] + [10]
new_data['text_color'] = ds.data['text_color'] + [RdYlBu3[i%3]]
new_data['text'] = ds.data['text'] + [str(i)]
ds.data = new_data

# create a callback that will add a number in a random location
def callback():
    global i, xs, ys, xe, ye

    # BEST PRACTICE --- update .data in one step with a new dict
    new_data = dict()
    new_data['x'] = ds.data['x'] + [random()*(xe-xs) + xs]
    new_data['y'] = ds.data['y'] + [random()*(ye-ys) + ys]
    new_data['text_color'] = ds.data['text_color'] + [RdYlBu3[i%3]]
    new_data['text'] = ds.data['text'] + [str(i)]
    ds.data = new_data

    i = i + 1

# create a callback that will add a number in a random location
def callback2():
    global i

    # BEST PRACTICE --- update .data in one step with a new dict
    new_data = {"x":[], "y":[], 'text_color':[], "text":[]}
    ds.data = new_data

    i = i + 1

def x_range_changed(*args):
    print(*args)


# add a button widget and configure with the call back
button = Button(label="Press Me")
button.on_click(callback)
# add a button widget and configure with the call back
button2 = Button(label="Press Me 2")
button2.on_click(callback2)

callback = CustomJS(args=dict(xr=p.x_range, yr=p.y_range, xs=xs, xe=xe, ys=ys, ye=ye), code="""
    function dist(a, b){
        return Math.abs( a - b );
    }
    if(dist(xs, xr.start) > 10 || dist(xe, xr.end) > 10 || dist(ys, yr.start) > 10 || dist(ye, yr.end) > 10)
    {
        s = document.location.href.split("?")[0] + "?xs=" + xr.start +
                                                "&xe=" + xr.end +
                                                "&ys=" + yr.start +
                                                "&ye=" + yr.end;
        //alert(s);
        document.location.href = s;
    }
""")
p.x_range.js_on_change('start', callback)
p.y_range.js_on_change('start', callback)

# put the button and plot in a layout and add to the document
curdoc().add_root(column(button, button2, p))