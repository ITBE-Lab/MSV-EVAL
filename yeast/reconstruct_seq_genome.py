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
    run_id = 3

    db_conn = DbConn({"SCHEMA": {"NAME": db_name}})
    call_table = SvCallTable(db_conn)
    jump_table = SvJumpTable(db_conn) # initialize jump table
    pack, _, _, _ = load_genomes(query_genome, reference_genome)
    seeds_list = call_table.calls_to_seeds(pack, run_id, True, 10)

    reconstructed_query_genome = call_table.reconstruct_sequenced_genome_from_seeds(seeds_list, pack)
    reconstructed_query_genome.store(reconstructed_query_genome_path + "/ma/genome")
    
    seeds_n_rects_reconstr = compute_seeds(MinimizerSeeding(ParameterSetManager()), reconstructed_query_genome_path,
                                           query_genome, 2, None)


    seeds_list_display = [(reconstructed_query_genome.start_of_sequence("reconstructed_" + name), seeds, []) for name, seeds, _ in seeds_list]

    seeds_n_rects = compute_seeds(MinimizerSeeding(ParameterSetManager()), query_genome,
                                  reference_genome, 2, filter_by_k_mer_set)

    out = []
    out.append(render_seeds_2(seeds_list_display, seeds_n_rects, query_genome, reference_genome))
    out.append(render_seeds_2(seeds_n_rects_reconstr, None, reconstructed_query_genome_path, query_genome))
    show(row(out))