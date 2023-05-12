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
    print("rendering", out_file)
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
                    elif dash == "dotted":
                        plot.circle(x=table[0][3],
                                    y=table[0][4],
                                    size=point_to_px(10),
                                    fill_color=color_scheme(color),
                                    legend_label=name,
                                    line_color=None,
                                    line_width=point_to_px(4))
                    elif dash == "dashed":
                        plot.square(x=table[0][3],
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

class NoneForNonExistant:
    def __init__(self, run_table):
        self.run_table = run_table
    
    def getId(self, name):
        if not self.run_table.hasName(name):
            return None
        return self.run_table.getId(name)

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
    dataset = "ufrj50816"
    prefix = "100nt-"
    run_table = NoneForNonExistant(SvCallerRunTable(DbConn({"SCHEMA": {"NAME": dataset}})))
    gt_lt_ten = run_table.getId("Ground Truth - Small <10 ")
    gt_ge_ten = run_table.getId("Ground Truth - Small >=10")
    gt_large = run_table.getId("Ground Truth - Large")
    recompute = True 
    run_ids = {
        gt_lt_ten: [run_table.getId("MA < 10nt 100nt"), 
                    run_table.getId("gridss-100 - Small < 10"), 
                    run_table.getId("manta-100 - Small < 10")],
        gt_ge_ten: [run_table.getId("MA [10,200)nt 100nt"), 
                    run_table.getId("gridss-100 - Small"), 
                    run_table.getId("manta-100 - Small")],
        gt_large:  [run_table.getId("MA >= 200nt 100nt"), 
                    run_table.getId("gridss-100 - Large"), 
                    run_table.getId("manta-100 - Large")],
    }
    if recompute:
        compute_accuracy_recall("ufrj50816", [0, 25], [gt_lt_ten, gt_ge_ten, gt_large], run_ids,
                                accuracy_recall_data_dir + "/" + prefix)

    render_multi_acc_recall_pic(
        [
            ("/" + prefix + dataset + "-" + str(run_ids[gt_lt_ten][0]) + "-0.tsv", "MSV 0nt", "blue", "solid"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_lt_ten][0]) + "-25.tsv", "MSV 25nt", "blue", "dotted"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_lt_ten][1]) + "-0.tsv", "Gridss 0nt", "orange", "solid"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_lt_ten][1]) + "-25.tsv", "Gridss 25nt", "orange", "dotted"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_lt_ten][2]) + "-0.tsv", "Manta 0nt", "green", "solid"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_lt_ten][2]) + "-25.tsv", "Manta 25nt", "green", "dotted"),
        ],
        "ufrj50816 [0,10) 100nt"
    )
    render_multi_acc_recall_pic(
        [
            ("/" + prefix + dataset + "-" + str(run_ids[gt_ge_ten][0]) + "-0.tsv", "MSV 0nt", "blue", "solid"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_ge_ten][0]) + "-25.tsv", "MSV 25nt", "blue", "dotted"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_ge_ten][1]) + "-0.tsv", "Gridss 0nt", "orange", "solid"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_ge_ten][1]) + "-25.tsv", "Gridss 25nt", "orange", "dotted"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_ge_ten][2]) + "-0.tsv", "Manta 0nt", "green", "solid"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_ge_ten][2]) + "-25.tsv", "Manta 25nt", "green", "dotted"),
        ],
        "ufrj50816 [10,200) 100nt"
    )
    render_multi_acc_recall_pic(
        [
            ("/" + prefix + dataset + "-" + str(run_ids[gt_large][0]) + "-0.tsv", "MSV 0nt", "blue", "solid"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_large][0]) + "-25.tsv", "MSV 25nt", "blue", "dotted"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_large][1]) + "-0.tsv", "Gridss 0nt", "orange", "solid"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_large][1]) + "-25.tsv", "Gridss 25nt", "orange", "dotted"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_large][2]) + "-0.tsv", "Manta 0nt", "green", "solid"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_large][2]) + "-25.tsv", "Manta 25nt", "green", "dotted"),
        ],
        "ufrj50816 [200,inf) 100nt"
    )

    prefix = "250nt-"
    run_ids = {
        gt_lt_ten: [run_table.getId("MA < 10nt 250nt"), 
                    run_table.getId("gridss-250 - Small < 10"), 
                    run_table.getId("manta-250 - Small < 10")],
        gt_ge_ten: [run_table.getId("MA [10,200)nt 250nt"), 
                    run_table.getId("gridss-250 - Small"), 
                    run_table.getId("manta-250 - Small")],
        gt_large:  [run_table.getId("MA >= 200nt 250nt"), 
                    run_table.getId("gridss-250 - Large"), 
                    run_table.getId("manta-250 - Large")],
    }
    if recompute:
        compute_accuracy_recall("ufrj50816", [0, 25], [gt_lt_ten, gt_ge_ten, gt_large], run_ids,
                                accuracy_recall_data_dir + "/" + prefix)

    render_multi_acc_recall_pic(
        [
            ("/" + prefix + dataset + "-" + str(run_ids[gt_lt_ten][0]) + "-0.tsv", "MSV 0nt", "blue", "solid"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_lt_ten][0]) + "-25.tsv", "MSV 25nt", "blue", "dotted"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_lt_ten][1]) + "-0.tsv", "Gridss 0nt", "orange", "solid"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_lt_ten][1]) + "-25.tsv", "Gridss 25nt", "orange", "dotted"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_lt_ten][2]) + "-0.tsv", "Manta 0nt", "green", "solid"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_lt_ten][2]) + "-25.tsv", "Manta 25nt", "green", "dotted"),
        ],
        "ufrj50816 [0,10) 250nt"
    )
    render_multi_acc_recall_pic(
        [
            ("/" + prefix + dataset + "-" + str(run_ids[gt_ge_ten][0]) + "-0.tsv", "MSV 0nt", "blue", "solid"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_ge_ten][0]) + "-25.tsv", "MSV 25nt", "blue", "dotted"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_ge_ten][1]) + "-0.tsv", "Gridss 0nt", "orange", "solid"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_ge_ten][1]) + "-25.tsv", "Gridss 25nt", "orange", "dotted"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_ge_ten][2]) + "-0.tsv", "Manta 0nt", "green", "solid"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_ge_ten][2]) + "-25.tsv", "Manta 25nt", "green", "dotted"),
        ],
        "ufrj50816 [10,200) 250nt"
    )
    render_multi_acc_recall_pic(
        [
            ("/" + prefix + dataset + "-" + str(run_ids[gt_large][0]) + "-0.tsv", "MSV 0nt", "blue", "solid"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_large][0]) + "-25.tsv", "MSV 25nt", "blue", "dotted"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_large][1]) + "-0.tsv", "Gridss 0nt", "orange", "solid"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_large][1]) + "-25.tsv", "Gridss 25nt", "orange", "dotted"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_large][2]) + "-0.tsv", "Manta 0nt", "green", "solid"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_large][2]) + "-25.tsv", "Manta 25nt", "green", "dotted"),
        ],
        "ufrj50816 [200,inf) 250nt"
    )

    prefix = "pacb-"
    run_ids = {
        gt_lt_ten: [run_table.getId("MA < 10nt PacBio")],
        gt_ge_ten: [run_table.getId("MA [10,200)nt PacBio")],
        gt_large:  [run_table.getId("MA >=200nt PacBio")],
    }
    if recompute:
        compute_accuracy_recall("ufrj50816", [0, 25], [gt_lt_ten, gt_ge_ten, gt_large], run_ids,
                                accuracy_recall_data_dir + "/" + prefix)
    render_multi_acc_recall_pic(
        [
            ("/" + prefix + dataset + "-" + str(run_ids[gt_lt_ten][0]) + "-0.tsv", "MSV 0nt", "blue", "solid"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_lt_ten][0]) + "-25.tsv", "MSV 25nt", "blue", "dotted"),
        ],
        "ufrj50816 [0,10) pacb"
    )
    render_multi_acc_recall_pic(
        [
            ("/" + prefix + dataset + "-" + str(run_ids[gt_ge_ten][0]) + "-0.tsv", "MSV 0nt", "blue", "solid"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_ge_ten][0]) + "-25.tsv", "MSV 25nt", "blue", "dotted"),
        ],
        "ufrj50816 [10,200) pacb"
    )
    render_multi_acc_recall_pic(
        [
            ("/" + prefix + dataset + "-" + str(run_ids[gt_large][0]) + "-0.tsv", "MSV 0nt", "blue", "solid"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_large][0]) + "-25.tsv", "MSV 25nt", "blue", "dotted"),
        ],
        "ufrj50816 [200,inf) pacb"
    )


    prefix = "realillumina-"
    run_ids = {
        gt_lt_ten: [run_table.getId("MA < 10nt realWorldIllumina"), 
                    run_table.getId("gridss-real-world - Small < 10"), 
                    run_table.getId("manta-real-world - Small < 10")],
        gt_ge_ten: [run_table.getId("MA [10,200)nt realWorldIllumina"), 
                    run_table.getId("gridss-real-world - Small"), 
                    run_table.getId("manta-real-world - Small")],
        gt_large:  [run_table.getId("MA >=200nt realWorldIllumina"), 
                    run_table.getId("gridss-real-world - Large"), 
                    run_table.getId("manta-real-world - Large")],
    }
    if recompute:
        run_ids = {
            gt_lt_ten: [#run_table.getId("gridss-real-world - Small < 10"), 
                        run_table.getId("manta-real-world - Small < 10")],
            gt_ge_ten: [#run_table.getId("gridss-real-world - Small"), 
                        run_table.getId("manta-real-world - Small")],
            gt_large:  [#run_table.getId("gridss-real-world - Large"), 
                        run_table.getId("manta-real-world - Large")],
        }
        #print(run_ids)
        compute_accuracy_recall("ufrj50816", [0, 25], [gt_lt_ten, gt_ge_ten, gt_large], run_ids,
                                accuracy_recall_data_dir + "/" + prefix)
    run_ids = {
        gt_lt_ten: [run_table.getId("MA < 10nt realWorldIllumina"), 
                    run_table.getId("gridss-real-world - Small < 10"), 
                    run_table.getId("manta-real-world - Small < 10")],
        gt_ge_ten: [run_table.getId("MA [10,200)nt realWorldIllumina"), 
                    run_table.getId("gridss-real-world - Small"), 
                    run_table.getId("manta-real-world - Small")],
        gt_large:  [run_table.getId("MA >=200nt realWorldIllumina"), 
                    run_table.getId("gridss-real-world - Large"), 
                    run_table.getId("manta-real-world - Large")],
    }
    render_multi_acc_recall_pic(
        [
            ("/" + prefix + dataset + "-" + str(run_ids[gt_lt_ten][0]) + "-0.tsv", "MSV 0nt", "blue", "solid"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_lt_ten][0]) + "-25.tsv", "MSV 25nt", "blue", "dotted"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_lt_ten][1]) + "-0.tsv", "Gridss 0nt", "orange", "solid"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_lt_ten][1]) + "-25.tsv", "Gridss 25nt", "orange", "dotted"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_lt_ten][2]) + "-0.tsv", "Manta 0nt", "green", "solid"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_lt_ten][2]) + "-25.tsv", "Manta 25nt", "green", "dotted"),
        ],
        "ufrj50816 [0,10) illumina real world"
    )
    render_multi_acc_recall_pic(
        [
            ("/" + prefix + dataset + "-" + str(run_ids[gt_ge_ten][0]) + "-0.tsv", "MSV 0nt", "blue", "solid"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_ge_ten][0]) + "-25.tsv", "MSV 25nt", "blue", "dotted"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_ge_ten][1]) + "-0.tsv", "Gridss 0nt", "orange", "solid"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_ge_ten][1]) + "-25.tsv", "Gridss 25nt", "orange", "dotted"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_ge_ten][2]) + "-0.tsv", "Manta 0nt", "green", "solid"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_ge_ten][2]) + "-25.tsv", "Manta 25nt", "green", "dotted"),
        ],
        "ufrj50816 [10,200) illumina real world"
    )
    render_multi_acc_recall_pic(
        [
            ("/" + prefix + dataset + "-" + str(run_ids[gt_large][0]) + "-0.tsv", "MSV 0nt", "blue", "solid"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_large][0]) + "-25.tsv", "MSV 25nt", "blue", "dotted"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_large][1]) + "-0.tsv", "Gridss 0nt", "orange", "solid"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_large][1]) + "-25.tsv", "Gridss 25nt", "orange", "dotted"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_large][2]) + "-0.tsv", "Manta 0nt", "green", "solid"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_large][2]) + "-25.tsv", "Manta 25nt", "green", "dotted"),
        ],
        "ufrj50816 [200,inf) illumina real world"
    )

    prefix = "realpacb-"
    run_ids = {
        gt_lt_ten: [
                        run_table.getId("MA < 10nt realWorldPacBio")
                    ],
        gt_ge_ten: [run_table.getId("MA [10,200)nt realWorldPacBio")],
        gt_large:  [run_table.getId("MA >=200nt realWorldPacBio")],
    }
    if recompute:
        compute_accuracy_recall("ufrj50816", [0, 25], [gt_lt_ten, gt_ge_ten, gt_large], run_ids,
                                accuracy_recall_data_dir + "/" + prefix)
    render_multi_acc_recall_pic(
        [
            ("/" + prefix + dataset + "-" + str(run_ids[gt_lt_ten][0]) + "-0.tsv", "MSV 0nt", "blue", "solid"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_lt_ten][0]) + "-25.tsv", "MSV 25nt", "blue", "dotted"),
        ],
        "ufrj50816 [0,10) pacb real world"
    )
    render_multi_acc_recall_pic(
        [
            ("/" + prefix + dataset + "-" + str(run_ids[gt_ge_ten][0]) + "-0.tsv", "MSV 0nt", "blue", "solid"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_ge_ten][0]) + "-25.tsv", "MSV 25nt", "blue", "dotted"),
        ],
        "ufrj50816 [10,200) pacb real world"
    )
    render_multi_acc_recall_pic(
        [
            ("/" + prefix + dataset + "-" + str(run_ids[gt_large][0]) + "-0.tsv", "MSV 0nt", "blue", "solid"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_large][0]) + "-25.tsv", "MSV 25nt", "blue", "dotted"),
        ],
        "ufrj50816 [200,inf) pacb real world"
    )

    prefix = "realoxfnano-"
    run_ids = {
        gt_lt_ten: [
                        run_table.getId("MA < 10nt oxfNano")
                    ],
        gt_ge_ten: [run_table.getId("MA [10,200)nt oxfNano")],
        gt_large:  [run_table.getId("MA >=200nt oxfNano")],
    }
    if recompute or True:
        compute_accuracy_recall("ufrj50816", [0, 25], [gt_lt_ten, gt_ge_ten, gt_large], run_ids,
                                accuracy_recall_data_dir + "/" + prefix)
    render_multi_acc_recall_pic(
        [
            ("/" + prefix + dataset + "-" + str(run_ids[gt_lt_ten][0]) + "-0.tsv", "MSV 0nt", "blue", "solid"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_lt_ten][0]) + "-25.tsv", "MSV 25nt", "blue", "dotted"),
        ],
        "ufrj50816 [0,10) oxfNano real world"
    )
    render_multi_acc_recall_pic(
        [
            ("/" + prefix + dataset + "-" + str(run_ids[gt_ge_ten][0]) + "-0.tsv", "MSV 0nt", "blue", "solid"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_ge_ten][0]) + "-25.tsv", "MSV 25nt", "blue", "dotted"),
        ],
        "ufrj50816 [10,200) oxfNano real world"
    )
    render_multi_acc_recall_pic(
        [
            ("/" + prefix + dataset + "-" + str(run_ids[gt_large][0]) + "-0.tsv", "MSV 0nt", "blue", "solid"),
            ("/" + prefix + dataset + "-" + str(run_ids[gt_large][0]) + "-25.tsv", "MSV 25nt", "blue", "dotted"),
        ],
        "ufrj50816 [200,inf) oxfNano real world"
    )
