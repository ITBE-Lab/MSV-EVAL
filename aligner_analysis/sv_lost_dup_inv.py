from aligner_analysis.binary_search_plot import *

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

def dup_inv(sv_size, gap_size, genome_section, ref_start, min_dist=50):
    points = []
    offset = random.randint(min_dist, (len(genome_section) - sv_size*2)-min_dist)
    g = str(genome_section)
    # before duplication
    read = g[:offset] 

    # section a (section to be duplicated)
    read += g[offset:offset + sv_size] 

    # section b (section that was duplicated)
    read += invert(g[offset:offset + sv_size])
    points.append((offset+sv_size, ref_start+offset+sv_size, True))
    points.append((offset+sv_size, ref_start+offset+sv_size, False))
    points.append((offset+sv_size*2, ref_start+offset, False))

    # after duplication
    read += g[offset + sv_size:len(genome_section)-sv_size]
    points.append((offset+sv_size*2, ref_start+offset+sv_size, True))

    return points, NucSeq(read)

def dup_inv_size(read_size, sv_size, gap_size):
    assert sv_size < read_size
    return read_size - sv_size

def main():
    params = ParameterSetManager()

    pack = Pack()
    pack.load(genome_dir + "/ma/genome")
    fm_index = FMIndex()
    fm_index.load(genome_dir + "/ma/genome")

    seeds_by_name, read_by_name = create_reads(pack, 1000, 100, lambda x,y: dup_inv(100, 20, x,y))
    path_sam = create_alignment(read_by_name, mm2, "mm2")
    print("Minimap 2 alignment:")
    comp = compare_alignment_from_file_paths(params, read_by_name, seeds_by_name, pack, path_sam)
    print("overlapped:", 100 * comp.nt_overlap / comp.nt_ground_truth, "% (nt)",
          100 * comp.amount_overlap / comp.amount_ground_truth, "% (seeds)")
    print("SMEMs:")
    comp = compare_seeds(params, read_by_name, seeds_by_name, fm_index, pack, True)
    print("overlapped:", 100 * comp.nt_overlap / comp.nt_ground_truth, "% (nt)",
          100 * comp.amount_overlap / comp.amount_ground_truth, "% (seeds)")

#main()
if False:
    binary_search_plot(dup_inv, "dup_inv_overlap")
    print_binary_search_plot_box_plot(file_name_in="dup_inv_overlap", title="Overlap - Duplication")

if True:
    accuracy_plot(dup_inv, dup_inv_size, "dup_inv_overlap")
    print_accuracy_plot(file_name_in="dup_inv_overlap", title="Accuracy - Duplication & Inversion")