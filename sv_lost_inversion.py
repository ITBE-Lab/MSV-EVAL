from binary_search_plot import *
import random

only_sv = True
no_sv = False

def invert(seq):
    ret = ""
    for c in seq[::-1]:
        ret += {
            "A": "T",
            "T": "A",
            "C": "G",
            "G": "C",
            "N": "N",

            "a": "t",
            "t": "a",
            "c": "g",
            "g": "c",
            "n": "n",
        }[c]
    return ret

def inversion(sv_size, gap_size, genome_section, ref_start, min_dist=50):
    seeds = Seeds()
    total_size = sv_size + gap_size
    offset = random.randint(min_dist, (len(genome_section) - total_size)-min_dist)
    total_size += offset
    g = str(genome_section)
    # before inversion
    if not only_sv:
        seeds.append(Seed(0, offset, ref_start, True))
    read = g[:offset] 

    # inverted secion
    if not no_sv:
        seeds.append(Seed(offset, sv_size, ref_start + total_size - 1, False))
    read += invert( g[offset + gap_size:total_size] )

    # after inversion
    if not only_sv:
        seeds.append(Seed(sv_size + offset, len(genome_section) - total_size, ref_start + total_size, True))
    read += g[total_size:len(genome_section)]

    return seeds, NucSeq(read)

def main():
    params = ParameterSetManager()

    pack = Pack()
    pack.load(genome_dir + "/ma/genome")
    fm_index = FMIndex()
    fm_index.load(genome_dir + "/ma/genome")

    seeds_by_name, read_by_name, gt_comp = create_reads(pack, 1000, 100, lambda x,y: inversion(500, 50, x, y))
    #path_sam = create_alignment(read_by_name, lambda x,y,z: mm2(x,y,z,True), "mm2")
    path_sam = create_alignment(read_by_name, mm2, "mm2")
    print("Minimap 2 alignment:")
    comp = compare_alignment_from_file_paths(params, read_by_name, seeds_by_name, pack, path_sam, gt_comp, True)
    print("overlapped:", 100 * comp.nt_overlap / comp.nt_ground_truth, "% (nt)",
          100 * comp.amount_overlap / comp.amount_ground_truth, "% (seeds)")
    print("SMEMs:")
    comp = compare_seeds(params, read_by_name, seeds_by_name, fm_index, pack, gt_comp, True)
    print("overlapped:", 100 * comp.nt_overlap / comp.nt_ground_truth, "% (nt)",
          100 * comp.amount_overlap / comp.amount_ground_truth, "% (seeds)")

#plot_quads(inversion)
#main()
if True:
    test_sets=[MM2TestSet("--splice"), SeedsTestSet(), NgmlrTestSet()]
    binary_search_plot(inversion, "inversion_overlap", test_sets=test_sets)
    #print_binary_search_plot("inversion_overlap", "SV Overlap - Inversion", test_sets=test_sets)