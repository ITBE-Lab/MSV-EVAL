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
    pack, _, _, ret_query_genome = load_genomes(query_genome, reference_genome)


    seeds_list = call_table.calls_to_seeds_by_id_auto(pack, run_id, True, 0)

    reconstructed_query_genome = call_table.reconstruct_sequenced_genome_from_seeds(seeds_list, pack)
    reconstructed_query_genome.store(reconstructed_query_genome_path + "/ma/genome")

    seeds_n_rects_reconstr = compute_seeds(reconstructed_query_genome_path, query_genome, db_name, 1)


    seeds_list_display = [(reconstructed_query_genome.start_of_sequence(name), seeds, []) for name, seeds, _ in seeds_list]

    seeds_n_rects = compute_seeds(query_genome, reference_genome, db_name, 1)

    out = []
    out.append(render_seeds_2(seeds_list_display, None, reconstructed_query_genome_path, reference_genome, title="reconstructed on reference"))
    out.append(render_seeds_2(seeds_n_rects, None, query_genome, reference_genome, title="assembly on reference"))
    out.append(render_seeds_2(seeds_n_rects_reconstr, None, reconstructed_query_genome_path, query_genome, title="reconstructed on assembly"))
    show(row(out))

    if True:
        print("name", "score", "matches", "missmatches", "indels", "indel ops", "% identity", sep="\t")
        for name, reconstr, (_, assembly) in zip(reconstructed_query_genome.contigNames(),
                                                reconstructed_query_genome.contigNucSeqs(),
                                                ret_query_genome):
            matches = 0
            mismatches = 0
            indels = 0
            indelops = 0
            l_total = max(len(reconstr), len(assembly))
            nw_alignment = runKsw(reconstr, assembly)
            for op, l in nw_alignment.data:
                if op == MatchType.match or op == MatchType.seed:
                    matches += l
                if op == MatchType.missmatch:
                    mismatches += l
                if op == MatchType.insertion or op == MatchType.deletion:
                    indels += l
                    indelops += 1
            iden = 100 * matches / l_total
            print(name, nw_alignment.get_score(), matches, mismatches, indels, indelops, iden, sep="\t")
