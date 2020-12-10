from aligner_analysis.binary_search_plot import *
import random


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
    points = []
    total_size = sv_size *2
    offset = random.randint(min_dist, (len(genome_section) - total_size)-min_dist)
    total_size += offset
    g = str(genome_section)

    # before inversion
    read = g[:offset] 
    points.append((offset, ref_start+offset, True))

    # inverted secion
    read += invert( g[offset + sv_size:total_size] )
    points.append((offset+sv_size, ref_start+offset+sv_size-1, False))
    points.append((offset, ref_start+total_size-1, False))

    # after inversion
    read += g[total_size:len(genome_section)]
    points.append((offset+sv_size, ref_start+total_size, True))

    return points, NucSeq(read)

def print_sam_file():
    params = ParameterSetManager()

    pack = Pack()
    pack.load(genome_dir + "/ma/genome")
    fm_index = FMIndex()
    fm_index.load(genome_dir + "/ma/genome")

    seeds_by_name, read_by_name, _ = create_reads(pack, 1000, 1, lambda x,y: inversion(100, 0, x, y))
    path_sam = create_alignment(read_by_name, lambda x,y,z: mm2(x,y,z,"-z 400,1"), "mm2-inversion-print_sam_file")

    sam_file = FileStreamFromPath(path_sam[0])
    file_reader = SamFileReader(params)
    lumper = SeedLumping(params)
    get_read_by_name = GetReadByName(params)
    align_to_seeds = AlignmentToSeeds(params)
    get_seeds_by_name = GetSeedsByReadName(params)
    idx = 0
    while not sam_file.eof():
        alignment = file_reader.execute(sam_file, pack, read_by_name)
        print(alignment.cigarString(pack))
        alignment_seeds = align_to_seeds.execute(alignment, pack)
        read = get_read_by_name.execute(alignment, read_by_name)
        lumped_data = lumper.execute( alignment_seeds, read, pack)
        printer = SeedPrinter(params, name_a="Alignment " + str(idx))
        gt_data = get_seeds_by_name.execute(read, seeds_by_name)
        printer.execute( alignment_seeds, gt_data )
        idx += 1

def main():
    params = ParameterSetManager()

    pack = Pack()
    pack.load(genome_dir + "/ma/genome")
    fm_index = FMIndex()
    fm_index.load(genome_dir + "/ma/genome")

    seeds_by_name, read_by_name = create_reads(pack, 1000, 1, lambda x,y: inversion(10, 0, x, y))

    if False:
        print("Minimap 2 alignment:")
        path_sam = create_alignment(read_by_name, lambda x,y,z: mm2(x,y,z,"-z 400,1"), "mm2-inversion-main")
        print(path_sam)
        comp = compare_alignment_from_file_paths(params, read_by_name, seeds_by_name, pack, path_sam, False)
        print("overlapped:", 100 * comp.nt_overlap / comp.nt_ground_truth, "% (nt)",
            100 * comp.amount_overlap / comp.amount_ground_truth, "% (seeds)")
    if False:
        print("MA alignment:")
        comp = MATestSet(render_one=False).test(params, seeds_by_name, read_by_name, fm_index,
                                               pack, "ma-inversion-main")
        print("overlapped:", 100 * comp.nt_overlap / comp.nt_ground_truth, "% (nt)",
            100 * comp.amount_overlap / comp.amount_ground_truth, "% (seeds)")
    print("reseeding:")
    compare_seeds(params, read_by_name, seeds_by_name, fm_index, pack, render_one=True)
    comp = compare_seeds(params, read_by_name, seeds_by_name, fm_index, pack)
    print("overlapped:", 100 * comp.nt_overlap / comp.nt_ground_truth, "% (nt)",
          100 * comp.amount_overlap / comp.amount_ground_truth, "% (seeds)")

#main()
#plot_quads(inversion)
#print_sam_file()

if True:
    accuracy_plot(inversion, filename_out="inversion_overlap")
    print_accuracy_plot(file_name_in="inversion_overlap", title="Accuracy - Inversion")