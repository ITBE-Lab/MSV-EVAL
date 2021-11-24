from bokeh.models import FactorRange
from bokeh.models import PrintfTickFormatter
from bokeh.plotting import figure, reset_output
from bokeh.transform import dodge
from bokeh.layouts import column
from bokeh.io import output_file, export_svgs, show, save
from sv_util.bokeh_style_helper import *
from sv_util.settings import *
from bokeh.models import NumeralTickFormatter
from MSV import *


def render_accuracy_n_recall_pic(tsv_file, outfile, y_range=(-0.051,1.051)):
    table = None
    num_gt = 0
    with open(tsv_file, "r") as in_file:
        lines = in_file.readlines()
        num_gt = float(lines[0][len("//|ground truth| = "):])
        table = [[float(y) for y in x.split("\t")] for x in lines[2:]]
    # supporting_reads, true_positives, num_entries, recall, accuracy = table[0]
    output_file(accuracy_recall_data_dir + "/bokeh_out_" + outfile + ".html")
    plot = figure(title=outfile, plot_height=300, plot_width=600, y_range=y_range)
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
    plot.output_backend = "svg"
    if show_plots:
        show(plot)
    if save_plots:
        save(plot)
    #export_svgs(plot, filename=accuracy_recall_data_dir + "/bokeh_out_" + outfile + ".svg")

def render_multi_acc_recall_pic(tsv_and_name_file, out_file, x_range=(-0.05,1.05),y_range=(-0.05,1.05)):
    output_file(accuracy_recall_data_dir + "/bokeh_out_" + out_file + ".html")
    plot = figure(title=out_file, plot_height=600, plot_width=600,x_range=x_range,y_range=y_range)
    plot.yaxis.axis_label = "accuracy [%]"
    plot.xaxis.axis_label = "recall [%]"
    plot.xaxis.formatter = NumeralTickFormatter(format='0%')
    plot.yaxis.formatter = NumeralTickFormatter(format='0%')
    for tsv_file, name, color, dash in tsv_and_name_file:
        try:
            with open(accuracy_recall_data_dir + tsv_file, "r") as in_file:
                lines = in_file.readlines()
                num_gt = float(lines[0][len("//|ground truth| = "):])
                table = [[float(y) for y in x.split("\t")] for x in lines[2:]]
                if len(table) == 1:
                    if dash == "solid":
                        plot.triangle(x=table[0][3],
                                    y=table[0][4],
                                    size=point_to_px(10),
                                    fill_color=color_scheme(color),
                                    legend_label=name,
                                    line_color=None,
                                    line_width=point_to_px(4))
                    else:
                        plot.circle(x=table[0][3],
                                    y=table[0][4],
                                    size=point_to_px(10),
                                    fill_color=color_scheme(color),
                                    legend_label=name,
                                    line_color=None,
                                    line_width=point_to_px(4))
                elif len(table) > 1:
                    plot.line(  x=[x[3] for x in table],
                                y=[x[4] for x in table],
                                legend_label=name,
                                line_color=color_scheme(color),
                                line_dash=dash,
                                line_width=point_to_px(4))
        except:
            pass
    style_plot(plot)
    plot.output_backend = "svg"
    if show_plots:
        show(plot)
    if save_plots:
        save(plot)


if __name__ == "__main__":
    #render_accuracy_n_recall_pic(accuracy_recall_data_dir + "/ufrj50816-1.tsv",
    #                             "ufrj50816-MA-Illumina-ge-10")
    #render_accuracy_n_recall_pic(accuracy_recall_data_dir + "/ufrj50816-2.tsv",
    #                             "ufrj50816-MA-PacBio")
    #render_accuracy_n_recall_pic(accuracy_recall_data_dir + "/ufrj50816-20.tsv",
    #                             "ufrj50816-MA-Illumina-st-10")
    #render_accuracy_n_recall_pic(accuracy_recall_data_dir + "/ufrj50816-3.tsv",
    #                             "ufrj50816-Gridss-PacBio")
    #render_accuracy_n_recall_pic(accuracy_recall_data_dir + "/ufrj50816-5.tsv",
    #                             "ufrj50816-Gridss-Illumina-ge-10")
    if True:
        gt_lt_ten = 36
        gt_ge_ten = 35
        gt_large = 33
        run_ids = {
            gt_lt_ten: [83],#[20, 32, 73, 65],
            gt_ge_ten: [82],#[1, 5, 72, 64],
            gt_large:  [],#[2, 3, 70, 62],
        }
        compute_accuracy_recall("ufrj50816", [0, 25], [gt_lt_ten, gt_ge_ten, gt_large], run_ids,
                                accuracy_recall_data_dir + "/")

    render_multi_acc_recall_pic(
        [
            ("/ufrj50816-20-0.tsv", "MSV 0nt", "blue", "solid"),
            ("/ufrj50816-20-25.tsv", "MSV 25nt", "blue", "dotted"),
            ("/ufrj50816-83-0.tsv", "MSV 0nt 100nt reads", "black", "solid"),
            ("/ufrj50816-83-25.tsv", "MSV 25nt 100nt reads", "black", "dotted"),
            ("/ufrj50816-32-0.tsv", "Gridss 0nt", "orange", "solid"),
            ("/ufrj50816-32-25.tsv", "Gridss 25nt", "orange", "dotted"),
            ("/ufrj50816-65-0.tsv", "Gridss 0nt 100nt reads", "red", "solid"),
            ("/ufrj50816-65-25.tsv", "Gridss 25nt 100nt reads", "red", "dotted"),
            ("/ufrj50816-81-0.tsv", "Manta 0nt", "green", "solid"),
            ("/ufrj50816-81-25.tsv", "Manta 25nt", "green", "dotted"),
        ],
        "ufrj50816 [0,10)"
    )
    render_multi_acc_recall_pic(
        [
            ("/ufrj50816-1-0.tsv", "MSV 0nt", "blue", "solid"),
            ("/ufrj50816-1-25.tsv", "MSV 25nt", "blue", "dotted"),
            ("/ufrj50816-82-0.tsv", "MSV 0nt 100nt reads", "black", "solid"),
            ("/ufrj50816-82-25.tsv", "MSV 25nt 100nt reads", "black", "dotted"),
            ("/ufrj50816-5-0.tsv", "Gridss 0nt", "orange", "solid"),
            ("/ufrj50816-5-25.tsv", "Gridss 25nt", "orange", "dotted"),
            ("/ufrj50816-64-0.tsv", "Gridss 0nt 100nt reads", "red", "solid"),
            ("/ufrj50816-64-25.tsv", "Gridss 25nt 100nt reads", "red", "dotted"),
            ("/ufrj50816-80-0.tsv", "Manta 0nt", "green", "solid"),
            ("/ufrj50816-80-25.tsv", "Manta 25nt", "green", "dotted"),
        ],
        "ufrj50816 [10,200)",
        x_range=(-0.05,0.55)
    )
    render_multi_acc_recall_pic(
        [
            ("/ufrj50816-2-0.tsv", "MSV 0nt", "blue", "solid"),
            ("/ufrj50816-2-25.tsv", "MSV 25nt", "blue", "dotted"),
            ("/ufrj50816-3-0.tsv", "Gridss 0nt", "orange", "solid"),
            ("/ufrj50816-3-25.tsv", "Gridss 25nt", "orange", "dotted"),
            ("/ufrj50816-62-0.tsv", "Gridss 0nt 100nt reads", "red", "solid"),
            ("/ufrj50816-62-25.tsv", "Gridss 25nt 100nt reads", "red", "dotted"),
            ("/ufrj50816-78-0.tsv", "Manta 0nt", "green", "solid"),
            ("/ufrj50816-78-25.tsv", "Manta 25nt", "green", "dotted"),
        ],
        "ufrj50816 [200,inf)"
    )

    #render_accuracy_n_recall_pic("/usr/home/markus/workspace/aligner/build.zeus/MSV/sv_visualization/ufrj50816-5.tsv",
    #                             "ufrj50816-PacBio-Simulated")
    #render_accuracy_n_recall_pic("/usr/home/markus/workspace/aligner/build.zeus/MSV/sv_visualization/ufrj50816-5.tsv",
    #                             "ufrj50816-PacBio-Simulated-zoom", x_range=(30, 60), y_range=(0.6,0.851))
    #render_accuracy_n_recall_pic("/usr/home/markus/workspace/aligner/build.zeus/MSV/sv_visualization/ufrj50816-6.tsv",
    #                             "ufrj50816-Illumina-Simulated")
    #render_accuracy_n_recall_pic("/usr/home/markus/workspace/aligner/build.zeus/MSV/sv_visualization/ufrj50816-6.tsv",
    #                             "ufrj50816-Illumina-Simulated-zoom", x_range=(10,70), y_range=(0.8,1))
    #render_accuracy_n_recall_pic(accuracy_recall_data_dir + "/ufrj50816-8.tsv",
    #                             "ufrj50816-Illumina-RealWorld")