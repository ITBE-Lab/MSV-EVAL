from svs_hidden_to_aligners.svs_lost_during_alignment import *
from sv_util.settings import *
from MA import *
from bokeh.models import FactorRange
from bokeh.models import PrintfTickFormatter
from bokeh.transform import dodge
from bokeh.layouts import column, gridplot
from bokeh.plotting import save
from bokeh.io import output_file, export_svgs
from sv_util.bokeh_style_helper import *
from bokeh.models import NumeralTickFormatter, FixedTicker
import subprocess

class MATestSet:
    def __init__(self, params=ParameterSetManager(), name="ma", render_one=False):
        params.by_name("Detect Small Inversions").set(True)
        self.params = params
        self._name = name
        self.render_one = render_one

    def test(self, params, seeds_by_name, read_by_name, fm_index, mm_index, pack, suffix, genome_section_by_name):
        path_sam = sv_hidden_to_aligners_data_dir + "/sam/" + self._name + suffix + ".sam"
        # write reads
        reads_path = sv_hidden_to_aligners_data_dir + "/reads/" + self._name + suffix + ".fasta"
        with open(reads_path, 'w') as fasta_file:
            for name, read in read_by_name:
                fasta_file.write(">" + name + "\n")
                fasta_file.write(str(read) + "\n")
        
        # align
        quick_align_paths([reads_path], human_genome_dir + "/ma/genome", self.params, path_sam)
        #print(reads_path)
        return compare_alignment_from_file_paths(params, read_by_name, seeds_by_name, pack, fm_index, [path_sam],
                                                 self.render_one)
    def name(self):
        return self._name
    def display_name(self):
        return "MA"

    def bokeh_func(self, plot, x, y, c, l):
        plot.x(x=x, y=y, line_color=c, legend_label=l,
                size=point_to_px(7), line_width=point_to_px(2))

    def color(self):
        return "orange"
    def color_light(self):
        return "yellow"

class MM2TestSet:
    def __init__(self, mm_extra="", name="mm2", presetting_override=None):
        self.mm2 = lambda x,y,z: mm2(x,y,z,extra=mm_extra, presetting_override=presetting_override)
        self._name = name
        self.c = "red"
        if len(mm_extra) > 0:
            self.c = "purple"
        elif not presetting_override is None:
            self.c = "turquise"

        self.c2 = "lightred"
        if len(mm_extra) > 0:
            self.c2 = "purple"
        elif not presetting_override is None:
            self.c = "turquise"

        self._display_name = "Minimap 2 Alignment"
        if len(mm_extra) > 0:
            self._display_name = "Minimap 2 Alignment extra"
    def test(self, params, seeds_by_name, read_by_name, fm_index, mm_index, pack, suffix, genome_section_by_name):
        path_sam = create_alignment(read_by_name, self.mm2, self._name + suffix)
        return compare_alignment_from_file_paths(params, read_by_name, seeds_by_name, pack, fm_index, path_sam)
    def name(self):
        return self._name
    def display_name(self):
        return self._display_name

    def bokeh_func(self, plot, x, y, c, l):
        plot.triangle(x=x, y=y, line_color=c, legend_label=l,
                size=point_to_px(7), line_width=point_to_px(2))

    def color(self):
        return self.c
    def color_light(self):
        return self.c2

class PBMM2TestSet:
    def __init__(self,name="pbmm2"):
        self.pbmm2 = lambda x,y,z: pbmm2(x,y,z)
        self._name = name
        self.c = "red"
        self.c2 = "lightred"
        self._display_name = "PB-Minimap 2 Alignment"

    def test(self, params, seeds_by_name, read_by_name, fm_index, mm_index, pack, suffix, genome_section_by_name):
        path_sam = create_alignment(read_by_name, self.pbmm2, self._name + suffix)
        return compare_alignment_from_file_paths(params, read_by_name, seeds_by_name, pack, fm_index, path_sam)
    def name(self):
        return self._name
    def display_name(self):
        return self._display_name

    def bokeh_func(self, plot, x, y, c, l):
        plot.triangle(x=x, y=y, line_color=c, legend_label=l,
                size=point_to_px(7), line_width=point_to_px(2))

    def color(self):
        return self.c
    def color_light(self):
        return self.c2

class NWTestSet:
    def __init__(self, name="SW"):
        self.sw = lambda x,y,z: sw(x,y,z)
        self._name = name
        self.c = "grey"#"purple"
        self.c2 = "pink"

        self._display_name = "NW Alignment"
    def test(self, params, seeds_by_name, read_by_name, fm_index, mm_index, pack, suffix, genome_section_by_name):
        return compare_nw(params, read_by_name, genome_section_by_name, seeds_by_name, pack)

    def name(self):
        return self._name
    def display_name(self):
        return self._display_name

    def bokeh_func(self, plot, x, y, c, l):
        plot.square(x=x, y=y, line_color=c, legend_label=l,
                    size=point_to_px(7), line_width=point_to_px(2))

    def color(self):
        return self.c

    def color_light(self):
        return self.c2

class NgmlrTestSet:
    def __init__(self):
        pass
    def test(self, params, seeds_by_name, read_by_name, fm_index, mm_index, pack, suffix, genome_section_by_name):
        path_sam = create_alignment(read_by_name, ngmlr, "ngmlr-" + suffix)
        return compare_alignment_from_file_paths(params, read_by_name, seeds_by_name, pack, fm_index, path_sam)
    def name(self):
        return "ngmlr"

    def display_name(self):
        return "Ngmlr Alignment"

    def bokeh_func(self, plot, x, y, c, l):
        plot.plus(x=x, y=y, line_color=c, legend_label=l,
                size=point_to_px(7), line_width=point_to_px(2))

    def color(self):
        return "green"
    def color_light(self):
        return "lightgreen"

class SeedsTestSet:
    def __init__(self, reseeding=False, mems=True):
        self.reseeding = reseeding
        self.mems = True

    def test(self, params, seeds_by_name, read_by_name, fm_index, mm_index, pack, suffix, genome_section_by_name):
        if self.mems:
            return compare_seeds(params, read_by_name, seeds_by_name, mm_index, pack, True, self.reseeding)
        else:
            return compare_seeds(params, read_by_name, seeds_by_name, fm_index, pack, False, self.reseeding)
    def name(self):
        if self.reseeding:
            return "remems"
        else:
            return "mems"
    def display_name(self):
        if self.reseeding:
            return "Reseeding"
        if self.mems:
            return "MEMs"
        return "Max. Spanning Seeding"

    def bokeh_func(self, plot, x, y, c, l):
        plot.circle(x=x, y=y, line_color=c, fill_color=None, legend_label=l,
                    size=point_to_px(7), line_width=point_to_px(2))

    def color(self):
        return "black"#"blue"
    def color_light(self):
        return "lightblue"

class GraphAligner:
    def __init__(self):
        pass

    def test(self, params, seeds_by_name, read_by_name, fm_index, mm_index, pack, suffix, genome_section_by_name):
        path_sam = create_alignment(read_by_name, graph_aligner, "graphaligner-" + suffix)
        return compare_alignment_from_file_paths(params, read_by_name, seeds_by_name, pack, fm_index, path_sam)

    def name(self):
        return "GraphAligner"

    def display_name(self):
        return "GraphAligner"

    def bokeh_func(self, plot, x, y, c, l):
        plot.circle(x=x, y=y, line_color=c, fill_color=None, legend_label=l,
                    size=point_to_px(7), line_width=point_to_px(2))

    def color(self):
        return "blue"
    def color_light(self):
        return "lightblue"

class GraphAligner2:
    def __init__(self):
        pass

    def test(self, params, seeds_by_name, read_by_name, fm_index, mm_index, pack, suffix, genome_section_by_name):
        path_sam = create_alignment(read_by_name, graph_aligner_2, "graphaligner50-" + suffix)
        return compare_alignment_from_file_paths(params, read_by_name, seeds_by_name, pack, fm_index, path_sam)

    def name(self):
        return "GraphAligner50"

    def display_name(self):
        return "GraphAligner50"

    def bokeh_func(self, plot, x, y, c, l):
        plot.circle(x=x, y=y, line_color=c, fill_color=None, legend_label=l,
                    size=point_to_px(7), line_width=point_to_px(2))

    def color(self):
        return "blue"
    def color_light(self):
        return "lightblue"

def print_n_write(s, f, fn):
    print(s, end="")
    f.write(s)

default_test_set = [
                        #PBMM2TestSet(),
                        NWTestSet(),
                        SeedsTestSet(),
                        MATestSet(),
                        MM2TestSet(),
                        NgmlrTestSet(),
                        GraphAligner(),
                        #MM2TestSet("-z 400,1 --splice -P", "mm2_extra"), 
                        #MM2TestSet(name="asm10", presetting_override="asm10"), 
                    ]
#default_test_set = [ MM2TestSet(), SeedsTestSet() ]
#default_test_set = [ GraphAligner2() ]


default_range = range(50, 500, 50)
#default_range = range(0, 250, 50)
#default_range = range(50, 500, 200)
#default_range = [450]

def accuracy_plot(sv_func, size_func=lambda x,y,z: x, filename_out="translocation_overlap",
                    test_sets=default_test_set, sv_sizes=default_range, read_sizes=[20000, 2000, 1000], 
                    do_disfigures=[True, False],
                    num_reads=1000
                    ):
    #return
    params = ParameterSetManager()
    params.set_selected("SV-PacBio")

    pack = Pack()
    pack.load(human_genome_dir + "/ma/genome")
    fm_index = FMIndex()
    fm_index.load(human_genome_dir + "/ma/genome")
    mm_index = MinimizerIndex(params, human_genome_dir + "/ma/genome.mmi")
    for read_size in read_sizes:
        for do_disfigure in do_disfigures:
            with open(sv_hidden_to_aligners_data_dir + "/" + str(read_size) + "-" 
                      + ("w_err" if do_disfigure else "err_free") + "-" + filename_out + ".tsv", "w") as file_out:
                # header of outfile
                print_n_write("read_size\ttest_set", file_out, filename_out)
                for sv_size in sv_sizes:
                    print_n_write("\t", file_out, filename_out)
                    print_n_write(str(sv_size), file_out, filename_out)
                print_n_write("\t# " + filename_out + (" disfigured" if do_disfigure else " error-free"), 
                              file_out, filename_out)
                print_n_write("\n", file_out, filename_out)

                # body of outfile
                read_cache = {}
                for test_set in test_sets:
                    print_n_write(str(read_size), file_out, filename_out)
                    print_n_write("\t", file_out, filename_out)
                    print_n_write(test_set.name(), file_out, filename_out)
                    for sv_size in sv_sizes:
                        print_n_write("\t", file_out, filename_out)
                        if not sv_func is None:
                            seeds_by_name, read_by_name, genome_section_by_name = \
                                        create_reads(pack, size_func(read_size, sv_size, 0),
                                                            num_reads,
                                                            lambda x,y: sv_func(sv_size, 0, x,y),
                                                    do_disfigure)
                        else:
                            seeds_by_name, read_by_name, genome_section_by_name = \
                                                create_scattered_read(pack, num_reads,
                                                                        0, sv_size, do_disfigure)
                        suffix = filename_out + "-sv_size=" + str(sv_size)
                        try:
                            acc = test_set.test(params, seeds_by_name, read_by_name, fm_index, mm_index, pack, suffix,
                                                genome_section_by_name)
                        except subprocess.TimeoutExpired:
                            acc = float("NaN")
                        print_n_write(str(acc), file_out, filename_out)
                    print_n_write("\n", file_out, filename_out)

def print_accuracy_plot(file_name_in="scattered_overlap", title="Overlap - Scattered read", 
                            test_sets_1=default_test_set, x_label="SV Size [nt]", save_svg=True,
                            read_sizes=[20000, 2000, 1000],
                            do_disfigures=[True, False],
                            ):
    #return
    for read_size in read_sizes:
        for do_disfigure in do_disfigures:
            with open(sv_hidden_to_aligners_data_dir + "/" + str(read_size)+"-"
                      + ("w_err" if do_disfigure else "err_free") + "-" + file_name_in + ".tsv", "r") as file_in:
                lines = file_in.readlines()
                header = lines[0]

                test_set_dict = {}
                for test_set in test_sets_1:
                    test_set_dict[test_set.name()] = test_set

                res = 4
                xs = [int(x) for x in lines[0].split("\t")[2:-1]]
                xs = [xs[0]] + xs + [xs[-1]]


                for line in lines[1:]:
                    cells = line.split("\t")
                    if cells[1] in test_set_dict:
                        test_set = test_set_dict[cells[1]]
                        output_file(sv_hidden_to_aligners_data_dir + "/bokeh_out_"+ str(read_size) + \
                                    ("w_err" if do_disfigure else "err_free") + "_" + file_name_in + "_" + test_set.name() + ".html")

                        ys = [0] + [float(x) for x in cells[2:]] + [0]

                        plot = figure(title=title + " " + test_set.display_name(), plot_height=600, plot_width=600,
                                      y_range=(0,1), x_range=(xs[0],xs[-1]))
                        plot.xaxis.axis_label = x_label
                        plot.yaxis.axis_label = "Recall [%]"
                        plot.yaxis.formatter = NumeralTickFormatter(format='0%')
                        plot.yaxis.bounds = (0,1)
                        plot.xaxis.bounds = (xs[0],xs[-1])
                        
                        plot.xaxis.ticker = FixedTicker(ticks=xs[1:-1])
                        plot.xgrid.ticker = FixedTicker(ticks=xs[1:-1])
                        yts = [0, .25, .5, .75, 1]
                        plot.yaxis.ticker = FixedTicker(ticks=yts)
                        plot.ygrid.ticker = FixedTicker(ticks=yts)


                        plot.patch(
                            x=xs,
                            y=ys,
                            fill_color=color_scheme(test_set.color()),
                            line_color=color_scheme(test_set.color()),
                        )

                        style_plot(plot)
                        plot.output_backend = "svg"
                        if show_plots:
                            show(plot)
                        if save_plots:
                            save(plot)
                            if save_svg:
                                export_svgs(plot, filename=sv_hidden_to_aligners_data_dir + "/bokeh_out_" + str(read_size) + "_" + test_set.name() + "_" + ("w_err" if do_disfigure else "err_free") + file_name_in + ".svg")