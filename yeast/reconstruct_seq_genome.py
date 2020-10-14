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
    run_id = 2

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

    if True:
        print("name", "score", "matches", "missmatches", "indels", "indel ops", "% identity", sep="\t")
        xs = []
        cx = [0]
        ys = []
        cy = [0]
        for x_start, name, reconstr, (y_start, assembly) in zip(
                                                reconstructed_query_genome.contigStarts(),
                                                reconstructed_query_genome.contigNames(),
                                                reconstructed_query_genome.contigNucSeqs(),
                                                ret_query_genome):
            x = 0
            y = 0
            xs.append(x_start)
            ys.append(y_start)
            cx.append(x_start + len(reconstr))
            cy.append(y_start + len(assembly))
            matches = 0
            mismatches = 0
            indels = 0
            indelops = 0
            l_total = max(len(reconstr), len(assembly))
            nw_alignment = runKsw(reconstr, assembly, 1000)
            for op, l in nw_alignment.data:
                if op == MatchType.match or op == MatchType.seed:
                    matches += l
                    x += l
                    y += l
                if op == MatchType.missmatch:
                    mismatches += l
                    x += l
                    y += l
                    # break line
                    xs.append(float("NaN"))
                    ys.append(float("NaN"))
                if op == MatchType.insertion or op == MatchType.deletion:
                    indels += l
                    indelops += 1
                    if op == MatchType.insertion:
                        x += l
                    else:
                        y += l
                xs.append(x + x_start)
                ys.append(y + y_start)
            xs.append(float("NaN"))
            ys.append(float("NaN"))
            iden = 100 * matches / l_total
            print(name, nw_alignment.get_score(), matches, mismatches, indels, indelops, iden, sep="\t")
        plot = figure(title="alignments reconstructed on assembly", plot_width=1000, plot_height=1000)
        for x in cx:
            plot.line(x=[x, x], y=[0, cy[-1]], color="black")
        for y in cy:
            plot.line(x=[0, cx[-1]], y=[y, y], color="black")
        plot.line(x=xs, y=ys)
        out.append(plot)

    show(row(out))
