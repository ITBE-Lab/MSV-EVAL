from svs_hidden_to_aligners.binary_search_plot import *

only_sv = True
no_sv = False
def translocation(sv_size, gap_size, genome_section, ref_start, min_dist=100):
    points = []
    total_size = sv_size*3
    offset = random.randint(min_dist, (len(genome_section) - total_size)-min_dist)
    total_size += offset
    g = str(genome_section)

    # before translocation
    read = g[:offset] 
    points.append((offset, ref_start+offset, True))

    # section a
    read += g[offset + sv_size *2:total_size] 
    points.append((offset+sv_size*2, ref_start+offset, True))
    points.append((offset+sv_size*3, ref_start+offset+sv_size, True))

    # within translocation
    read += g[offset + sv_size:offset + sv_size *2]
    points.append((offset+sv_size, ref_start+offset+sv_size, True))
    points.append((offset+sv_size*2, ref_start+offset+sv_size*2, True))

    # section b
    read += g[offset:offset + sv_size]
    points.append((offset, ref_start+offset+sv_size*2, True))
    points.append((offset+sv_size, ref_start+offset+sv_size*3, True))

    # after tranlocation
    read += g[total_size:]
    points.append((total_size, ref_start+total_size, True))

    return points, NucSeq(read)

def main():
    params = ParameterSetManager()

    pack = Pack()
    pack.load(genome_dir + "/ma/genome")
    fm_index = FMIndex()
    fm_index.load(genome_dir + "/ma/genome")
    mm_index = MinimizerIndex(params, genome_dir + "/ma/genome.mmi")

    seeds_by_name, read_by_name = create_reads(pack, 1000, 1, lambda x,y: translocation(10, 0, x,y))
    if False:
        path_sam = create_alignment(read_by_name, mm2, "mm2")
        print("Minimap 2 alignment:")
        comp = compare_alignment_from_file_paths(params, read_by_name, seeds_by_name, pack, path_sam)
        print("overlapped:", 100 * comp.nt_overlap / comp.nt_ground_truth, "% (nt)",
            100 * comp.amount_overlap / comp.amount_ground_truth, "% (seeds)")
    print("reseeding:")
    compare_seeds(params, read_by_name, seeds_by_name, mm_index, pack, render_one=True)
    comp = compare_seeds(params, read_by_name, seeds_by_name, mm_index, pack)
    print("overlapped:", 100 * comp.nt_overlap / comp.nt_ground_truth, "% (nt)",
          100 * comp.amount_overlap / comp.amount_ground_truth, "% (seeds)")

#main()

if True:
    accuracy_plot(translocation, filename_out="translocation_overlap")
    print_accuracy_plot(file_name_in="translocation_overlap", title="Accuracy - Translocation")