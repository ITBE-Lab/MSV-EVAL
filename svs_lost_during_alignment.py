from MA import *
import random
from os_aligners import *
from bokeh.plotting import figure, show
from bokeh_style_helper import *
from bokeh.plotting import ColumnDataSource
from bokeh.layouts import column, row, grid
from bokeh.models.tools import HoverTool

genome_dir = "/MAdata/genome/human/GRCh38.p12-chr1"
data_dir = "/MAdata/sv_lost_during_alignment"

def choice_adj_size(l, total_len):
    x = random.randrange(total_len)
    idx = 0
    while x >= len(l[idx][0]):
        x -= len(l[idx][0])
        idx += 1
    return l[idx]

def from_to(contig, f, t):
    s = ""
    for x in range(f, f + t):
        s += contig[x]
    return NucSeq(s)

def create_reads(pack, size, amount, func_get_seeds_and_read):
    read_by_name = ReadByName()
    seeds_by_name = SeedsByName()
    contigs = [(x, y) for x, y in zip(pack.contigSeqs(), pack.contigStarts()) if len(x) > size]
    total_len = sum(len(x) for x, _ in contigs)
    def read_and_seeds():
        contig, contig_start = choice_adj_size(contigs, total_len)
        start = random.randrange(len(contig) - size)
        genome_section = from_to(contig, start, size)
        return func_get_seeds_and_read(genome_section, start + contig_start)
    for idx in range(amount):
        seeds, read = read_and_seeds()
        if 'n' in str(read) or 'N' in str(read):
            continue
        read.name = "read" + str(idx)
        read_by_name.append(read)
        seeds_by_name.append(seeds, read.name)
    return seeds_by_name, read_by_name


def compare(params, ground_truth, data, reads, pack_pledge, unlock_targets=None):
    collect = CollectSeedSetComps(params)
    compare_module = CompareSeedSets(params)
    lumper = SeedLumping(params)

    res = VectorPledge()
    for idx in range(params.get_num_threads()):
        lumped_g_t = promise_me(lumper, ground_truth[idx], reads[idx], pack_pledge)
        lumped_data = promise_me(lumper, data[idx], reads[idx], pack_pledge)
        comp = promise_me(compare_module, lumped_g_t, lumped_data)
        empty = promise_me(collect, comp)
        if unlock_targets is None:
            unlock = empty
        else:
            unlock = promise_me(UnLock(params, unlock_targets[idx]), empty)
        res.append(unlock)

    res.simultaneous_get(params.get_num_threads())
    return collect.cpp_module.collection


def compare_seeds(params, reads_by_name, seeds_by_name, fm_index, pack):
    splitter = NucSeqSplitter(params)
    lock = Lock(params)
    reads_by_name_pledge = Pledge()
    reads_by_name_pledge.set(reads_by_name)
    pack_pledge = Pledge()
    pack_pledge.set(pack)
    fm_index_pledge = Pledge()
    fm_index_pledge.set(fm_index)
    seeds_by_name_pledge = Pledge()
    seeds_by_name_pledge.set(seeds_by_name)

    seeding_module = BinarySeeding(params)
    extract_seeds = ExtractSeeds(params)
    get_seeds_by_name = GetSeedsByReadName(params)

    reads_vec = ContainerVectorNucSeq()
    for name, read in reads_by_name:
        reads_vec.append(read)
    reads_vec_pledge = Pledge()
    reads_vec_pledge.set(reads_vec)

    ground_truth = []
    data = []
    reads = []
    unlock_targets = []
    read = promise_me(splitter, reads_vec_pledge)
    for _ in range(params.get_num_threads()):
        locked_read = promise_me(lock, read)
        unlock_targets.append(locked_read)
        reads.append(locked_read)
        segments = promise_me(seeding_module, fm_index_pledge, locked_read)
        seeds = promise_me(extract_seeds, segments, fm_index_pledge, locked_read)
        data.append(seeds)
        ground_truth_seeds = promise_me(get_seeds_by_name, locked_read, seeds_by_name_pledge)
        ground_truth.append(ground_truth_seeds)

    return compare(params, ground_truth, data, reads, pack_pledge, unlock_targets)

def compare_alignment(params, reads_by_name_pledge, seeds_by_name, alignments, pack_pledge, unlock_targets=None):
    align_to_seeds = AlignmentToSeeds(params)
    get_seeds_by_name = GetSeedsByName(params)
    get_read_by_name = GetReadByName(params)
    seeds_by_name_pledge = Pledge()
    seeds_by_name_pledge.set(seeds_by_name)
    ground_truth = []
    data = []
    reads = []
    for idx in range(params.get_num_threads()):
        alignment_seeds = promise_me(align_to_seeds, alignments[idx])
        data.append(alignment_seeds)
        ground_truth_seeds = promise_me(get_seeds_by_name, alignments[idx], seeds_by_name_pledge)
        ground_truth.append(ground_truth_seeds)
        read = promise_me(get_read_by_name, alignments[idx], reads_by_name_pledge)
        reads.append(read)

    return compare(params, ground_truth, data, reads, pack_pledge, unlock_targets)



def compare_alignment_from_file_queue(params, reads_by_name, seeds_by_name, pack, queue_pledge):
    file_reader = SamFileReader(params)
    queue_picker = FilePicker(params)
    queue_placer = FileAlignmentPlacer(params)
    lock = Lock(params)
    reads_by_name_pledge = Pledge()
    reads_by_name_pledge.set(reads_by_name)
    pack_pledge = Pledge()
    pack_pledge.set(pack)

    alignments = []
    locked_files = []
    for _ in range(params.get_num_threads()):
        picked_file = promise_me(queue_picker, queue_pledge)
        locked_file = promise_me(lock, picked_file)
        alignment_ = promise_me(file_reader, locked_file, pack_pledge, reads_by_name_pledge)
        alignment = promise_me(queue_placer, alignment_, locked_file, queue_pledge)
        locked_files.append(locked_file)
        alignments.append(alignment)

    return compare_alignment(params, reads_by_name_pledge, seeds_by_name, alignments, pack_pledge, locked_files)


def compare_alignment_from_file_paths(params, reads_by_name, seeds_by_name, pack, file_paths):
    if file_paths is None:
        return None
    file_queue = FileQueue()
    for string in file_paths:
        file_queue.add(FileStreamFromPath(string))
    queue_pledge = Pledge()
    queue_pledge.set(file_queue)
    return compare_alignment_from_file_queue(params, reads_by_name, seeds_by_name, pack, queue_pledge)

def create_alignment(read_by_name, aligner, sam_name):
    reads_path = data_dir + "/reads/" + sam_name + ".fasta"
    with open(reads_path, 'w') as fasta_file:
        for name, read in read_by_name:
            fasta_file.write(">" + name + "\n")
            fasta_file.write(str(read) + "\n")

    json_dict = {"reference_path":genome_dir}
    read_json = {"technology":"pb", "name":"n/a", "fasta_file":reads_path}
    path_sam = data_dir + "/mm_sam/" + sam_name + ".sam"
    aligner(read_json, path_sam, json_dict)

    return [path_sam]

def no_sv(genome_section, ref_start):
    seeds = Seeds()
    seeds.append(Seed(0, len(genome_section), ref_start, True))
    return seeds, genome_section

only_sv = True
no_sv = False
def translocation(sv_size, gap_size, genome_section, ref_start):
    seeds = Seeds()
    total_size = sv_size*2 + gap_size
    offset = (len(genome_section) - total_size)//2
    total_size += offset
    g = str(genome_section)
    # before translocation
    if not only_sv:
        seeds.append(Seed(0, offset, ref_start, True))
    read = g[:offset] 

    # section a
    if not no_sv:
        seeds.append(Seed(offset + sv_size + gap_size, sv_size, ref_start + offset, True))
    read += g[offset + sv_size + gap_size:total_size] 

    # within translocation
    if not only_sv:
        seeds.append(Seed(offset + sv_size, gap_size, ref_start + offset + sv_size, True))
    read += g[offset + sv_size:offset + sv_size + gap_size]

    # section b
    if not no_sv:
        seeds.append(Seed(offset, sv_size, ref_start + offset + sv_size + gap_size, True))
    read += g[offset:offset + sv_size]

    # after tranlocation
    if not only_sv:
        seeds.append(Seed(total_size, len(genome_section) - total_size, ref_start + total_size, True))
    read += g[total_size:]

    return seeds, NucSeq(read)

def main():
    params = ParameterSetManager()

    pack = Pack()
    pack.load(genome_dir + "/ma/genome")
    fm_index = FMIndex()
    fm_index.load(genome_dir + "/ma/genome")

    seeds_by_name, read_by_name = create_reads(pack, 1000, 100, lambda x,y: translocation(300, 10, x,y))
    path_sam = create_alignment(read_by_name, mm2, "mm2")
    print("Minimap 2 alignment:")
    comp = compare_alignment_from_file_paths(params, read_by_name, seeds_by_name, pack, path_sam)
    print("overlapped:", 100 * comp.nt_overlap / comp.nt_ground_truth, "% (nt)",
          100 * comp.amount_overlap / comp.amount_ground_truth, "% (seeds)")
    print("SMEMs:")
    comp = compare_seeds(params, read_by_name, seeds_by_name, fm_index, pack)
    print("overlapped:", 100 * comp.nt_overlap / comp.nt_ground_truth, "% (nt)",
          100 * comp.amount_overlap / comp.amount_ground_truth, "% (seeds)")

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

def binary_search_plot(filename_out="translocation_overlap", gap_size_range=range(0, 81, 2),
                       overlap_percentages=[0.05, 0.45, 0.95], test_sets=[
                           MM2TestSet(), SeedsTestSet()
                       ], sv_size_range=(0, 200)):
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
                        seeds_by_name, read_by_name = create_reads(pack, 1000, 1000,
                                                            lambda x,y: translocation(sv_size, gap_size, x,y))
                        read_cache[(sv_size, gap_size)] = (seeds_by_name, read_by_name)
                    suffix = "sv_size=" + str(sv_size) + "-gap_size=" + str(gap_size)
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
                    
                    sv_size = bin_search(sv_size_range[0], sv_size_range[1], o_p, search_val_cache, func=test)
                    if sv_size == sv_size_range[1]:
                        sv_size = None
                    print_n_write(str(sv_size), file_out)
                    print_n_write("\t", file_out)
                print_n_write("\n", file_out)

def print_binary_search_plot(file_name_in="translocation_overlap", test_sets=[MM2TestSet(), SeedsTestSet()],
                            sv_size_range=(180, 0), x_range=(-5, 85)):
    
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

        plot = figure(title="SV Overlap - Translocation", y_range=sv_size_range, x_range=x_range)
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
            plot.line(x=[0,1], y=[sv_size_range[1]-10, sv_size_range[1]-10], line_width=point_to_px(4), 
                      legend_label=color_strings + " overlap", color="black", line_dash=dash)
        for test_set in test_sets:
            plot.line(x=[0,1], y=[sv_size_range[1]-10, sv_size_range[1]-10], line_width=point_to_px(4), 
                      legend_label=test_set.display_name(), color=color_scheme(test_set.color()))
        plot.legend.location = "bottom_left"
        show(plot)

binary_search_plot()
print_binary_search_plot()