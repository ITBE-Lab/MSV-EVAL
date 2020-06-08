from aligner_analysis.binary_search_plot import *
import random


def deletion(sv_size, gap_size, genome_section, ref_start, min_dist=50):
    points = []
    total_size = sv_size
    offset = random.randint(min_dist, (len(genome_section) - total_size)-min_dist)
    g = str(genome_section)

    # before deletion
    read = g[:offset] 
    points.append((offset, ref_start+offset, True))

    # after deletion
    read += g[offset+sv_size:]
    points.append((offset, ref_start+offset+sv_size, True))

    return points, NucSeq(read)

def del_size(read_size, sv_size, gap_size):
    return read_size + sv_size


def insertion(sv_size, gap_size, genome_section, ref_start, min_dist=50):
    points = []
    offset = random.randint(min_dist, (len(genome_section))-min_dist)
    g = str(genome_section)

    # before insertion
    read = g[:offset] 
    points.append((offset, ref_start+offset, True))

    #insertion
    for _ in range(sv_size):
        read += random.choice(["a","c","g","t"])

    # after insertion
    read += g[offset:]
    points.append((offset+sv_size, ref_start+offset, True))

    return points, NucSeq(read)

def ins_size(read_size, sv_size, gap_size):
    assert sv_size < read_size
    return read_size - sv_size


def main():
    params = ParameterSetManager()

    pack = Pack()
    pack.load(genome_dir + "/ma/genome")
    fm_index = FMIndex()
    fm_index.load(genome_dir + "/ma/genome")
    mm_index = MinimizerIndex(params, genome_dir + "/ma/genome.mmi")

    for sv_size in range(10, 800, 100):
        print("sv_size:", sv_size)
        seeds_by_name, read_by_name = create_reads(pack, 1000-sv_size, 100, lambda x,y: insertion(sv_size, 0, x, y))

        if True:
            path_sam = create_alignment(read_by_name, lambda x,y,z: mm2(x,y,z), "mm2-indel-main")
            print(path_sam)
            comp = compare_alignment_from_file_paths(params, read_by_name, seeds_by_name, pack, fm_index, path_sam,
                                                    False)
            print("Minimap 2 alignment:", 100 * comp, "%")
        if True:
            comp = MATestSet(render_one=False).test(params, seeds_by_name, read_by_name, fm_index, mm_index,
                                                pack, "ma-indel-main")
            print("MA alignment:", 100 * comp, "%")

#main()

if True:
    accuracy_plot(deletion, del_size, "deletion_overlap")
    print_accuracy_plot(file_name_in="deletion_overlap", title="Accuracy - Deletion")

    accuracy_plot(insertion, ins_size, "insertion_overlap")
    print_accuracy_plot(file_name_in="insertion_overlap", title="Accuracy - Insertion")