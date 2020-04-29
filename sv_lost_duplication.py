from binary_search_plot import *

only_sv = True
no_sv = False
def duplication(sv_size, gap_size, genome_section, ref_start, min_dist=50):
    seeds = Seeds()
    total_size = sv_size*2 + gap_size
    offset = random.randint(min_dist, (len(genome_section) - total_size)-min_dist)
    total_size += offset
    g = str(genome_section)
    # before duplication
    if not only_sv:
        seeds.append(Seed(0, offset, ref_start, True))
    read = g[:offset] 

    # section a
    if not only_sv:
        seeds.append(Seed(offset, sv_size, ref_start + offset, True))
    read += g[offset:offset + sv_size] 

    # between duplication
    if not only_sv:
        seeds.append(Seed(offset + sv_size, gap_size, ref_start + offset + sv_size, True))
    read += g[offset + sv_size:offset + sv_size + gap_size]

    # section b
    if not no_sv:
        seeds.append(Seed(offset + sv_size + gap_size, sv_size, ref_start + offset, True))
    read += g[offset:offset + sv_size]

    # after duplication
    if not only_sv:
        seeds.append(Seed(total_size, len(genome_section) - total_size, ref_start + offset + sv_size + gap_size, True))
    read += g[offset + sv_size + gap_size:len(genome_section)-sv_size]

    return seeds, NucSeq(read)

def main():
    params = ParameterSetManager()

    pack = Pack()
    pack.load(genome_dir + "/ma/genome")
    fm_index = FMIndex()
    fm_index.load(genome_dir + "/ma/genome")

    seeds_by_name, read_by_name = create_reads(pack, 1000, 100, lambda x,y: duplication(100, 20, x,y))
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
if True:
    #binary_search_plot(duplication, "duplication_overlap", sv_size_max=1000, read_size=5000, gap_size_range=range(1,500,25))
    print_binary_search_plot("duplication_overlap", "SV Overlap - Duplication", sv_size_max=1000)