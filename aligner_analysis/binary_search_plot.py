from aligner_analysis.svs_lost_during_alignment import *
from MA import *
from bokeh.models import FactorRange
from bokeh.models import PrintfTickFormatter
from bokeh.transform import dodge
from bokeh.layouts import column
from bokeh.io import output_file, export_svgs
from sv_util.bokeh_style_helper import *
class MATestSet:
    def __init__(self, params=ParameterSetManager(), name="ma", render_one=False):
        params.by_name("Detect Small Inversions").set(True)
        self.params = params
        self._name = name
        self.render_one = render_one

    def test(self, params, seeds_by_name, read_by_name, fm_index, mm_index, pack, suffix, genome_section_by_name):
        path_sam = data_dir + "/sam/" + self._name + suffix + ".sam"
        # write reads
        reads_path = data_dir + "/reads/" + self._name + suffix + ".fasta"
        with open(reads_path, 'w') as fasta_file:
            for name, read in read_by_name:
                fasta_file.write(">" + name + "\n")
                fasta_file.write(str(read) + "\n")
        # align
        quick_align_paths([reads_path], genome_dir + "/ma/genome", self.params, path_sam)
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
    def __init__(self, mm_extra="", name="mm2"):
        self.mm2 = lambda x,y,z: mm2(x,y,z,extra=mm_extra)
        self._name = name
        self.c = "red"
        if len(mm_extra) > 0:
            self.c = "black"

        self.c2 = "lightred"
        if len(mm_extra) > 0:
            self.c2 = "grey"

        self._display_name = "Minimap 2 Alignment"
        if len(mm_extra) > 0:
            self._display_name = "Minimap 2 Alignment extra sensitive"
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

class SwTestSet:
    def __init__(self, name="SW"):
        self.sw = lambda x,y,z: sw(x,y,z)
        self._name = name
        self.c = "purple"
        self.c2 = "pink"

        self._display_name = "SW Alignment"
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
    def __init__(self, reseeding, mems=True):
        self.reseeding = reseeding
        self.mems = True

    def test(self, params, seeds_by_name, read_by_name, fm_index, mm_index, pack, suffix, genome_section_by_name):
        if self.mems:
            return compare_seeds(params, read_by_name, seeds_by_name, mm_index, pack, True, self.reseeding)
        else:
            return compare_seeds(params, read_by_name, seeds_by_name, fm_index, pack, False, self.reseeding)
    def name(self):
        return "seeds"
    def display_name(self):
        if self.reseeding:
            return "Reseeding"
        if self.mems:
            return "MEMS"
        return "Max. Spanning Seeding"

    def bokeh_func(self, plot, x, y, c, l):
        plot.circle(x=x, y=y, line_color=c, fill_color=None, legend_label=l,
                    size=point_to_px(7), line_width=point_to_px(2))

    def color(self):
        return "blue"
    def color_light(self):
        return "lightblue"

def print_n_write(s, f):
    print(s, end="")
    f.write(s)

default_test_set = [SwTestSet(), SeedsTestSet(False, True), MATestSet(), MM2TestSet(),
                    MM2TestSet("-z 400,1 --splice -P", "extra_sensitive"), NgmlrTestSet()]
#default_test_set = [SwTestSet()]

def accuracy_plot(sv_func, size_func=lambda x,y,z: x, filename_out="translocation_overlap",
                    test_sets=default_test_set, sv_sizes=range(25, 501, 25), read_size=2000, num_reads=1000):
    params = ParameterSetManager()
    params.set_selected("SV-PacBio")

    pack = Pack()
    pack.load(genome_dir + "/ma/genome")
    fm_index = FMIndex()
    fm_index.load(genome_dir + "/ma/genome")
    mm_index = MinimizerIndex(params, genome_dir + "/ma/genome.mmi")

    with open(data_dir + "/" + filename_out + ".tsv", "w") as file_out:
        # header of outfile
        print_n_write("gap_size\ttest_set", file_out)
        for sv_size in sv_sizes:
            print_n_write("\t", file_out)
            print_n_write(str(sv_size), file_out)
        print_n_write("\n", file_out)

        # body of outfile
        read_cache = {}
        for test_set in test_sets:
            print_n_write("0", file_out)
            print_n_write("\t", file_out)
            print_n_write(test_set.name(), file_out)
            for sv_size in sv_sizes:
                print_n_write("\t", file_out)
                if not sv_func is None:
                    seeds_by_name, read_by_name, genome_section_by_name = \
                                create_reads(pack, size_func(read_size, sv_size, 0),
                                                    num_reads,
                                                    lambda x,y: sv_func(sv_size, 0, x,y))
                else:
                    seeds_by_name, read_by_name, genome_section_by_name = \
                                        create_scattered_read(pack, num_reads,
                                                                0, sv_size)
                suffix = filename_out + "-sv_size=" + str(sv_size)
                acc = test_set.test(params, seeds_by_name, read_by_name, fm_index, mm_index, pack, suffix,
                                    genome_section_by_name)
                print_n_write(str(acc), file_out)
            print_n_write("\n", file_out)

def print_accuracy_plot(file_name_in="scattered_overlap", title="Overlap - Scattered read", 
                            test_sets=default_test_set, x_label="SV Size", save_svg=True):
    with open(data_dir + "/" + file_name_in + ".tsv", "r") as file_in:
        output_file(data_dir + "/bokeh_out_" + file_name_in + ".html")
        lines = file_in.readlines()
        header = lines[0]

        test_set_dict = {}
        for test_set in test_sets:
            test_set_dict[test_set.name()] = test_set

        xs = [int(x) for x in lines[0].split("\t")[2:]]

        plot = figure(title=title, plot_width=800, y_range=(-0.05, 1.05), x_range=(-10,510))
        plot.xaxis.axis_label = x_label
        plot.yaxis.axis_label = "Accuracy"
        funcs = [
            plot.x,
            plot.circle,
            plot.cross,
            plot.square,
            plot.triangle,
            plot.diamond,
            plot.circle_x,
            plot.circle_cross,
        ]
        for line, func in zip(lines[1:],funcs):
            cells = line.split("\t")
            test_set = test_set_dict[cells[1]]
            plot.line(x=xs, y=[float(x) for x in cells[2:]], line_color=color_scheme(test_set.color()),
                      line_width=point_to_px(4),
                      legend_label=test_set.display_name())
            func(x=xs, y=[float(x) for x in cells[2:]], line_color=color_scheme(test_set.color()),
                 line_width=point_to_px(2), fill_color=None,
                 size=point_to_px(8), legend_label=test_set.display_name())
        plot.legend.location = "bottom_right"
        style_plot(plot)
        show(plot)
        if save_svg:
            plot.output_backend = "svg"
            export_svgs(plot, filename=data_dir + "/bokeh_out_" + file_name_in + ".svg")