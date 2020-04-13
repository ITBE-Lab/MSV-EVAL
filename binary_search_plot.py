from svs_lost_during_alignment import *


def plot_gap_size():
    params = ParameterSetManager()

    pack = Pack()
    pack.load(genome_dir + "/ma/genome")
    fm_index = FMIndex()
    fm_index.load(genome_dir + "/ma/genome")

    print("gap_size","mm2","max. sp", sep="\t")
    for gap_size in range(0, 50, 5):
        seeds_by_name, read_by_name = create_reads(pack, 1000, 100, lambda x,y: translocation(100, gap_size, x,y))
        path_sam = create_alignment(read_by_name, mm2, "mm2-gap_size=" + str(gap_size))
        comp = compare_alignment_from_file_paths(params, read_by_name, seeds_by_name, pack, path_sam)
        comp2 = compare_seeds(params, read_by_name, seeds_by_name, fm_index, pack)
        print(gap_size, 100 * comp.nt_overlap / comp.nt_ground_truth, 100 * comp2.nt_overlap / comp2.nt_ground_truth, sep="\t")

def plot_quads():
    params = ParameterSetManager()

    pack = Pack()
    pack.load(genome_dir + "/ma/genome")
    fm_index = FMIndex()
    fm_index.load(genome_dir + "/ma/genome")

    rects_align = {
        "x": [],
        "y": [],
        "fill_color": [],
        "val": [],
    }
    rects_seeds = {
        "x": [],
        "y": [],
        "fill_color": [],
        "val": [],
    }
    w = 5
    h = 5
    for sv_size in range(3, 200, h):
        print("sv_size", sv_size)
        for gap_size in range(3, 70, w):
            seeds_by_name, read_by_name = create_reads(pack, 1000, 100,
                                                        lambda x,y: translocation(sv_size, gap_size, x,y))
            path_sam = create_alignment(read_by_name, mm2, "mm2-sv_size=" + str(sv_size) + "-gap_size=" + str(gap_size))
            comp = compare_alignment_from_file_paths(params, read_by_name, seeds_by_name, pack, path_sam)
            rects_align["y"].append(sv_size)
            rects_align["x"].append(gap_size)
            c = format(light_spec_approximation(comp.nt_overlap / comp.nt_ground_truth))
            rects_align["fill_color"].append(c)
            rects_align["val"].append(100 * comp.nt_overlap / comp.nt_ground_truth)
            comp2 = compare_seeds(params, read_by_name, seeds_by_name, fm_index, pack)
            rects_seeds["y"].append(sv_size)
            rects_seeds["x"].append(gap_size)
            c2 = format(light_spec_approximation(comp2.nt_overlap / comp2.nt_ground_truth))
            rects_seeds["fill_color"].append(c2)
            rects_seeds["val"].append(100 * comp2.nt_overlap / comp2.nt_ground_truth)

    plot1 = figure(title="Minimap 2 Alignment", y_range=(200,0))
    plot1.rect(x="x", y="y", width=w, height=h, color="fill_color", source=ColumnDataSource(rects_align))
    plot1.add_tools(HoverTool(tooltips=[("overlap", "@val%"),]))
    plot1.yaxis.axis_label = "SV Size"
    plot1.xaxis.axis_label = "Gap Size"

    plot2 = figure(title="Max. Spanning Seeds", y_range=(200,0))
    plot2.rect(x="x", y="y", width=w, height=h, color="fill_color", source=ColumnDataSource(rects_seeds))
    plot2.add_tools(HoverTool(tooltips=[("overlap", "@val%"),]))
    plot2.yaxis.axis_label = "SV Size"
    plot2.xaxis.axis_label = "Gap Size"

    show(column(plot1, plot2))

class MM2TestSet:
    def __init__(self):
        pass
    def test(self, params, seeds_by_name, read_by_name, fm_index, pack, suffix):
        path_sam = create_alignment(read_by_name, mm2, "mm2-" + suffix)
        comp = compare_alignment_from_file_paths(params, read_by_name, seeds_by_name, pack, path_sam)
        return comp.nt_overlap / comp.nt_ground_truth
    def name(self):
        return "mm2"
    def display_name(self):
        return "Minimap 2 Alignment"

    def bokeh_func(self, plot, x, y, c, l):
        plot.x(x=x, y=y, line_color=c, legend_label=l,
                size=point_to_px(7), line_width=point_to_px(2))

    def color(self):
        return "red"

class NgmlrTestSet:
    def __init__(self):
        pass
    def test(self, params, seeds_by_name, read_by_name, fm_index, pack, suffix):
        path_sam = create_alignment(read_by_name, ngmlr, "ngmlr-" + suffix)
        comp = compare_alignment_from_file_paths(params, read_by_name, seeds_by_name, pack, path_sam)
        return comp.nt_overlap / comp.nt_ground_truth
    def name(self):
        return "ngmlr"

    def display_name(self):
        return "Ngmlr Alignment"

    def bokeh_func(self, plot, x, y, c, l):
        plot.x(x=x, y=y, line_color=c, legend_label=l,
                size=point_to_px(7), line_width=point_to_px(2))

    def color(self):
        return "green"

class SeedsTestSet:
    def __init__(self):
        pass
    def test(self, params, seeds_by_name, read_by_name, fm_index, pack, suffix):
        comp = compare_seeds(params, read_by_name, seeds_by_name, fm_index, pack)
        return comp.nt_overlap / comp.nt_ground_truth
    def name(self):
        return "seeds"
    def display_name(self):
        return "Max. Spanning Seeding"

    def bokeh_func(self, plot, x, y, c, l):
        plot.circle(x=x, y=y, line_color=c, fill_color=None, legend_label=l,
                    size=point_to_px(7), line_width=point_to_px(2))

    def color(self):
        return "blue"

def print_n_write(s, f):
    print(s, end="")
    f.write(s)

def binary_search_plot(sv_func, filename_out="translocation_overlap", gap_size_range=range(0, 81, 2),
                       overlap_percentages=[0.05, 0.45, 0.95], test_sets=[
                           MM2TestSet(), SeedsTestSet()
                       ], sv_size_max=200, read_size=1000):
    params = ParameterSetManager()

    pack = Pack()
    pack.load(genome_dir + "/ma/genome")
    fm_index = FMIndex()
    fm_index.load(genome_dir + "/ma/genome")

    with open(data_dir + "/" + filename_out + ".tsv", "w") as file_out:
        # header of outfile
        print_n_write("gap_size\ttest_set\t", file_out)
        for o_p in overlap_percentages:
            print_n_write(str(100*o_p), file_out)
            print_n_write("%\t", file_out)
        print_n_write("\n", file_out)

        # body of outfile
        for gap_size in gap_size_range:
            read_cache = {}
            for test_set in test_sets:
                print_n_write(str(gap_size), file_out)
                print_n_write("\t", file_out)
                print_n_write(test_set.name(), file_out)
                print_n_write("\t", file_out)
                def test(sv_size):
                    if (sv_size, gap_size) in read_cache:
                        seeds_by_name, read_by_name = read_cache[(sv_size, gap_size)]
                    else:
                        seeds_by_name, read_by_name = create_reads(pack, read_size, 1000,
                                                            lambda x,y: sv_func(sv_size, gap_size, x,y))
                        read_cache[(sv_size, gap_size)] = (seeds_by_name, read_by_name)
                    suffix = "-" + filename_out + "-sv_size=" + str(sv_size) + "-gap_size=" + str(gap_size)
                    return test_set.test(params, seeds_by_name, read_by_name, fm_index, pack, suffix)

                search_val_cache = {}
                for o_p in overlap_percentages:
                    def bin_search(bottom, top, search_val, cache, func):
                        if bottom == top:
                            return top
                        x = (top + bottom) // 2
                        if x in cache:
                            val = cache[x]
                        else:
                            val = func(x)
                            cache[x] = val
                        if search_val < val:
                            return bin_search(bottom, x, search_val, cache, func)
                        if search_val > val:
                            if x == bottom:
                                return top
                            return bin_search(x, top, search_val, cache, func)
                        return x
                    
                    sv_size = bin_search(0, sv_size_max, o_p, search_val_cache, func=test)
                    if sv_size == sv_size_max:
                        sv_size = None
                    print_n_write(str(sv_size), file_out)
                    print_n_write("\t", file_out)
                print_n_write("\n", file_out)

def print_binary_search_plot(file_name_in="translocation_overlap", title="SV Overlap - Translocation", 
                            test_sets=[MM2TestSet(), SeedsTestSet()],
                            sv_size_max=180, x_range=(-5, 85)):
    
    with open(data_dir + "/" + file_name_in + ".tsv", "r") as file_in:
        lines = file_in.readlines()
        header = lines[0]

        # this ignores the first two columns (gap_size and test_set)
        color_strings = header.split("\t")[2:-1]
        dashes = [(5,20), (20,5), (20,0)]
        #colors = [percentage_scale(x) for x in color_strings]
        #colors = [format(light_spec_approximation(float(x[:-1])/100)) for x in color_strings]

        data = {}
        for line in lines[1:]:
            cells = line.split("\t")
            gap_size = cells[0]
            test_set = cells[1]
            y_values = cells[2:-1]
            if not test_set in data:
                data[test_set] = {"x":[],"ys":{}}
            data[test_set]["x"].append(gap_size)
            for idx, y_value in enumerate(y_values):
                if not idx in data[test_set]["ys"]:
                    data[test_set]["ys"][idx] = []
                data[test_set]["ys"][idx].append(y_value)

        test_set_dict = {}
        for test_set in test_sets:
            test_set_dict[test_set.name()] = test_set

        plot = figure(title=title, y_range=(sv_size_max, 0), x_range=x_range)
        plot.yaxis.axis_label = "SV Size"
        plot.xaxis.axis_label = "Gap Size"
        for test_set_name, test_set_data in data.items():
            test_set = test_set_dict[test_set_name]
            color = color_scheme(test_set.color())
            x = test_set_data["x"]
            for idx, y in test_set_data["ys"].items():
                dash = dashes[idx]
                plot.line(x=x, y=y, line_color=color, line_width=point_to_px(4), line_dash=dash)
                #test_set.bokeh_func(plot, x, y, color, test_set.display_name())

        for color_strings, dash in zip(color_strings, dashes):
            plot.line(x=[0,1], y=[-10, -10], line_width=point_to_px(4), 
                      legend_label=color_strings + " overlap", color="black", line_dash=dash)
        for test_set in test_sets:
            plot.line(x=[0,1], y=[-10, -10], line_width=point_to_px(4), 
                      legend_label=test_set.display_name(), color=color_scheme(test_set.color()))
        plot.legend.location = "bottom_left"
        show(plot)