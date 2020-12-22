from bokeh.models import FactorRange
from bokeh.models import PrintfTickFormatter
from bokeh.plotting import figure, reset_output
from bokeh.transform import dodge
from bokeh.layouts import column
from bokeh.io import output_file, export_svgs, show
from sv_util.bokeh_style_helper import *
from bokeh.models import NumeralTickFormatter


data_dir = "/MAdata/sv_caller_analysis/accuracy_recall"

def render_accuracy_n_recall_pic(tsv_file, outfile, x_range=(-11,271), y_range=(-0.051,1.051)):
    table = None
    num_gt = 0
    with open(tsv_file, "r") as in_file:
        lines = in_file.readlines()
        num_gt = float(lines[0][len("//|ground truth| = "):])
        table = [[float(y) for y in x.split("\t")] for x in lines[2:]]
    # supporting_reads, true_positives, num_entries, recall, accuracy = table[0]
    output_file(data_dir + "/bokeh_out_" + outfile + ".html")
    plot = figure(title=outfile, plot_height=300, plot_width=600, y_range=y_range, x_range=x_range)
    plot.yaxis.formatter = NumeralTickFormatter(format='0%')
    plot.xaxis.axis_label = "supporting reads"
    symbol_every = 10
    plot.line(  x=[x[0] for x in table],
                y=[x[3] for x in table],
                legend_label="recall",
                line_color=color_scheme("blue"),
                line_width=point_to_px(4))
    plot.x(     x=[x[0] for x in table][::symbol_every] + [table[-1][0]],
                y=[x[3] for x in table][::symbol_every] + [table[-1][3]],
                size=13,
                legend_label="recall",
                line_color=color_scheme("blue"),
                line_width=point_to_px(2))
    plot.line(  x=[x[0] for x in table],
                y=[x[4] for x in table],
                legend_label="accuracy",
                line_color=color_scheme("green"),
                line_width=point_to_px(4))
    plot.cross( x=[x[0] for x in table][::symbol_every] + [table[-1][0]],
                y=[x[4] for x in table][::symbol_every] + [table[-1][4]],
                size=13,
                legend_label="accuracy",
                line_color=color_scheme("green"),
                fill_color=None,
                line_width=point_to_px(2))
    style_plot(plot)
    show(plot)
    plot.output_backend = "svg"
    export_svgs(plot, filename=data_dir + "/bokeh_out_" + outfile + ".svg")


if __name__ == "__main__":
    #render_accuracy_n_recall_pic("/usr/home/markus/workspace/aligner/build.zeus/MSV/sv_visualization/ufrj50816-5.tsv",
    #                             "ufrj50816-PacBio-Simulated")
    #render_accuracy_n_recall_pic("/usr/home/markus/workspace/aligner/build.zeus/MSV/sv_visualization/ufrj50816-5.tsv",
    #                             "ufrj50816-PacBio-Simulated-zoom", x_range=(30, 60), y_range=(0.6,0.851))
    #render_accuracy_n_recall_pic("/usr/home/markus/workspace/aligner/build.zeus/MSV/sv_visualization/ufrj50816-6.tsv",
    #                             "ufrj50816-Illumina-Simulated")
    #render_accuracy_n_recall_pic("/usr/home/markus/workspace/aligner/build.zeus/MSV/sv_visualization/ufrj50816-6.tsv",
    #                             "ufrj50816-Illumina-Simulated-zoom", x_range=(10,70), y_range=(0.8,1))
    render_accuracy_n_recall_pic("/usr/home/markus/workspace/aligner/build.zeus/MSV/sv_visualization/ufrj50816-7.tsv",
                                 "ufrj50816-PacBio-RealWorld")
    render_accuracy_n_recall_pic("/usr/home/markus/workspace/aligner/build.zeus/MSV/sv_visualization/ufrj50816-8.tsv",
                                 "ufrj50816-Illumina-RealWorld")