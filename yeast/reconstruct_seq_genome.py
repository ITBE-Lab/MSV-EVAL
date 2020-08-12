from yeast.display_mems import *


#query_genome = "knowlesiStrain"
#query_genome = "UFRJ50816-chrVII-section"
query_genome = genome_dir + "UFRJ50816"
reconstructed_query_genome_path = "/MAdata/genome/reconstructed/yeast/UFRJ50816"
#query_genome = "YPS138-chrVII-section"
#reference_genome = "YPS138-chrVII-section"
reference_genome = genome_dir + "YPS138"
#reference_genome = "vivax"
if __name__ == "__main__":
    db_name = "UFRJ50816_test_reconstruct"
    run_id = 1

    db_conn = DbConn({"SCHEMA": {"NAME": db_name}})
    call_table = SvCallTable(db_conn)
    jump_table = SvJumpTable(db_conn) # initialize jump table
    pack, _, _, _ = load_genomes(query_genome, reference_genome)

    #extract_contigs_from = [ # without order
    #    ("chrI", True),
    #    ("chrII", True),
    #    ("chrIII", True),
    #    ("chrIV", True),
    #    ("chrXI", True),
    #    ("chrVI", True),
    #    ("chrVII", True),
    #    ("chrVIII", True),
    #    ("chrXV", True),
    #    ("chrX", True),
    #    ("chrXIV", True),
    #    ("chrXIV", False),
    #    ("chrXII", True),
    #    ("chrIX", False),
    #    ("chrXVI", True),
    #]

    extract_contigs_from = [
        ("chrI", True, "rec_chrI"),
        ("chrII", True, "rec_chrII"),
        ("chrIII", True, "rec_chrIII"),
        ("chrXI", True, "rec_chrIV"),
        ("chrVIII", True, "rec_chrV"),
        ("chrVI", True, "rec_chrVI"),
        ("chrVII", True, "rec_chrVII"),
        ("chrVIII", True, "rec_chrVIII"),
        ("chrIII", False, "rec_chrIX"),
        ("chrX", True, "rec_chrX"),
        ("chrIV", True, "rec_chrXI"),
        ("chrXIV", True, "rec_chrXII"),
        ("chrXIV", False, "rec_chrXIII"),
        ("chrXII", True, "rec_chrXIV"),
        ("chrXV", True, "rec_chrXV"),
        ("chrXIII", False, "rec_chrXVI"),
    ]
    seeds_list = call_table.calls_to_seeds_by_id(pack, run_id, True, 10, extract_contigs_from)

    reconstructed_query_genome = call_table.reconstruct_sequenced_genome_from_seeds(seeds_list, pack)
    reconstructed_query_genome.store(reconstructed_query_genome_path + "/ma/genome")
    
    seeds_n_rects_reconstr = compute_seeds(MinimizerSeeding(ParameterSetManager()), reconstructed_query_genome_path,
                                           query_genome, 2, None)


    seeds_list_display = [(reconstructed_query_genome.start_of_sequence(name), seeds, []) for (_, seeds, _), (_, f, name) in zip(seeds_list, extract_contigs_from)]

    seeds_n_rects = compute_seeds(MinimizerSeeding(ParameterSetManager()), query_genome,
                                  reference_genome, 2, filter_by_k_mer_set)

    out = []
    out.append(render_seeds_2(seeds_list_display, None, reconstructed_query_genome_path, reference_genome, title="reconstructed on reference"))
    out.append(render_seeds_2(seeds_n_rects, None, query_genome, reference_genome, title="assembly on reference"))
    out.append(render_seeds_2(seeds_n_rects_reconstr, None, reconstructed_query_genome_path, query_genome, title="reconstructed on assembly"))
    show(row(out))