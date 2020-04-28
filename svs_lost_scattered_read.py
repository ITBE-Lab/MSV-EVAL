from binary_search_plot import *












def plot_scatter(filename_out="scattered_overlap",  test_sets=[MM2TestSet(), SeedsTestSet()],
                 sv_size_range=range(1, 200, 10),
                     read_size=1000, num_reads=1000, num_scatters=10):
    params = ParameterSetManager()

    pack = Pack()
    pack.load(genome_dir + "/ma/genome")
    fm_index = FMIndex()
    fm_index.load(genome_dir + "/ma/genome")

    with open(data_dir + "/" + filename_out + ".tsv", "w") as file_out:
        # header of outfile
        print_n_write("sv_size\ttest_set\t", file_out)
        for idx in range(num_scatters+1):
            print_n_write("%found " + str(idx) + "\t", file_out)
        print_n_write("\n", file_out)

        # body of outfile
        for sv_size in sv_size_range:
            seeds_by_name, read_by_name, gt_comp = create_scattered_read(pack, read_size, num_reads,
                                                                         num_scatters, sv_size)
            for test_set in test_sets:
                print_n_write(str(sv_size), file_out)
                print_n_write("\t", file_out)
                print_n_write(test_set.name(), file_out)
                print_n_write("\t", file_out)
                
                suffix = filename_out + "-sv_size=" + str(sv_size)
                comp = test_set.test(params, seeds_by_name, read_by_name, fm_index, pack, suffix, gt_comp)
                for idx in range(num_scatters+1):
                    if idx in comp.seeds_found:
                        print_n_write(str(100*comp.seeds_found[idx]/num_reads) + "\t", file_out)
                    else:
                        print_n_write("0\t", file_out)
                print_n_write("\n", file_out)


##
# @brief computes the median
# @details
# find the index where the difference of the sum of all elements below and all elements above are smaller than the
# element itself
def q2(l):
    below = 0
    above = sum(l[1:])
    idx = 0
    while abs(below - above) > l[idx]:
        below += l[idx]
        idx += 1
        above -= l[idx]
    return idx

def draw_scatter_plot(filename_out="scattered_overlap", test_sets=[MM2TestSet(), SeedsTestSet()]):
    data = {}
    for test_set in test_sets:
        data[test_set.name()] = {"x": [], "min": [], "q2": [], "max": []}
    with open(data_dir + "/" + filename_out + ".tsv", "r") as file_in:
        body = file_in.readlines()[1:]
        for line in body:
            cells = line.split("\t")[:-1]
            percentages = [float(x) for x in cells[2:]]
            name = cells[1]
            data[name]["x"].append(float(cells[0]))
            data[name]["min"].append(min(idx if val > 0 else len(percentages)-1 for idx, val in enumerate(percentages)))
            data[name]["q2"].append(q2(percentages))
            data[name]["max"].append(max(idx if val > 0 else 0 for idx, val in enumerate(percentages)))
    
    test_set_dict = {}
    for test_set in test_sets:
        test_set_dict[test_set.name()] = test_set
    plot = figure(title="Overlap - Scattered read", x_axis_type="log")
    plot.xaxis.axis_label = "SV Size"
    plot.yaxis.axis_label = "Num Scatters"
    for test_set_name, test_set_data in data.items():
        test_set = test_set_dict[test_set_name]
        color = color_scheme(test_set.color())
        plot.patch(x=test_set_data["x"] + test_set_data["x"][::-1], 
                   y=test_set_data["min"] + test_set_data["max"][::-1], 
                   fill_color=color_scheme(test_set.color_light()), line_color=None)
        plot.line(x=test_set_data["x"], y=test_set_data["q2"], line_color=color, line_width=point_to_px(4),
                  legend_label=test_set_name)
    plot.legend.location = "top_left"
    show(plot)


def main():
    params = ParameterSetManager()

    pack = Pack()
    pack.load(genome_dir + "/ma/genome")
    fm_index = FMIndex()
    fm_index.load(genome_dir + "/ma/genome")

    num_scatters = 10
    amount = 100
    sv_size = 1000
    print_one = False

    seeds_by_name, read_by_name = create_scattered_read(pack, 1000, amount, num_scatters, sv_size)
    #path_sam = create_alignment(read_by_name, lambda x,y,z: mm2(x,y,z,True), "mm2")
    path_sam = create_alignment(read_by_name, mm2, "mm2")
    print("Minimap 2 alignment:")
    comp = compare_alignment_from_file_paths(params, read_by_name, seeds_by_name, pack, path_sam, print_one)
    print("overlapped:", 100 * comp.nt_overlap / comp.nt_ground_truth, "% (nt)",
          100 * comp.amount_overlap / comp.amount_ground_truth, "% (seeds)")
    for key in comp.seeds_found:
        val = comp.seeds_found[key] # pybind11 std::map type
        print("found", key, "of", num_scatters, 100*val/amount, "% of times")
    print("SMEMs:")
    comp = compare_seeds(params, read_by_name, seeds_by_name, fm_index, pack, print_one)
    print("overlapped:", 100 * comp.nt_overlap / comp.nt_ground_truth, "% (nt)",
          100 * comp.amount_overlap / comp.amount_ground_truth, "% (seeds)")
    for key in comp.seeds_found:
        val = comp.seeds_found[key] # pybind11 std::map type
        print("found", key, "of", num_scatters, 100*val/amount, "% of times")

mm_extra = ""
mm_extra += "--splice" # Enable the splice alignment mode.
#mm_extra += " -P"      # Retain all chains 

test_sets=[MM2TestSet(mm_extra), SeedsTestSet(), NgmlrTestSet()]
#plot_scatter(num_reads=100, sv_size_range=[2**x for x in range(14)], test_sets=test_sets, num_scatters=100)
draw_scatter_plot(test_sets=test_sets)
