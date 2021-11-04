from bokeh.events import DoubleTap
from bokeh.models.callbacks import CustomJS

def point_to_px(val):
    return val * 1.33333


def color_scheme(idx):
    colors = {
        "blue": "#4e75a3",
        "lightblue": "#b3c5da",
        "green": "#7d9263",
        "lightgreen" : "#c9d3be",
        "red" : "#c00000",
        "lightred" : "#ffafaf",
        "orange" : "#dd7e0e",
        "yellow" : "#dfcf04",
        "pink" : "#b55475",
        "purple" : "#7153a1",
        "black": "#000000",
        "grey": "#555555",
        "teal": "#008080",
        "lightteal": "##00ffff",
    }
    return colors[idx]

def percentage_scale(s):
    colors= {
        "5.0%": "yellow",
        "45.0%": "orange",
        "95.0%": "red",
    }
    return color_scheme(colors[s])

def light_spec_approximation(x):
    #map input [0, 1] to wavelength [350, 645]
    w = 370 + x * (645-370)
    r = 0.0
    g = 0.0
    b = 0.0
    if w < 440:
        r = -(w - 440.) / (440. - 380.)
        b = 1.0
    elif w >= 440 and w < 490:
        g = (w - 440.) / (490. - 440.)
        b = 1.0
    elif w >= 490 and w < 510:
        g = 1.0
        b = -(w - 510.) / (510. - 490.)
    elif w >= 510 and w < 580:
        r = (w - 510.) / (580. - 510.)
        g = 1.0
    elif w >= 580 and w < 645:
        r = 1.0
        g = -(w - 645.) / (645. - 580.)
    elif w >= 645:
        r = 1.0

    #intensity
    i = 1.0
    if w > 650:
        i = .3 + .7*(780-w)/(780-650)
    elif w < 420:
        i = .3 + .7*(w-380)/(420-380)

    #gamma
    m = .8

    return (i*r**m, i*g**m, i*b**m)

def format(rgb):
        def clamp(x):
            return max(0, min(x, 255))
        red, green, blue = rgb
        return "#{0:02x}{1:02x}{2:02x}".format(clamp(int(red * 255)), clamp(int(green * 255)),
                                               clamp(int(blue * 255)))

def style_plot(plot):
    plot.grid.grid_line_width = point_to_px(2)
    plot.outline_line_width = point_to_px(2)
    plot.axis.axis_line_width = point_to_px(2)
    plot.axis.major_tick_line_width = point_to_px(2)
    plot.axis.minor_tick_line_width = point_to_px(1)
    plot.axis.axis_label_text_color = "black"
    plot.axis.axis_label_text_font_size = "20pt"
    plot.axis.axis_label_text_font = "calibri"
    plot.axis.axis_label_text_font_style = "normal"
    plot.axis.major_label_text_color = "black"
    plot.axis.major_label_text_font_size = "20pt"
    plot.axis.major_label_text_font = "calibri"
    plot.axis.major_tick_out = 10
    plot.axis.major_tick_in = 4
    plot.axis.minor_tick_out = 7
    plot.legend.label_text_font = "calibri"
    plot.legend.label_text_font_size = "20pt"
    plot.legend.label_text_color = "black"
    plot.legend.label_height = 0
    plot.legend.label_text_line_height = 0
    plot.legend.spacing = 0
    plot.legend.background_fill_alpha = 1
    plot.legend.border_line_alpha = 1
    plot.legend.border_line_color = "black"
    plot.legend.border_line_width = point_to_px(1)
    plot.title.text_font_size  = "20pt"
    plot.title.text_color = "black"
    plot.title.text_font = "calibri"
    plot.js_on_event(DoubleTap, CustomJS(args=dict(other=plot.legend[0]),
             code="other.visible = !other.visible;"
    ))