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
    db_name = "UFRJ50816"
    run_ids = [1,2,3]

    db_conn = DbConn({"SCHEMA": {"NAME": db_name}})
    call_table = SvCallTable(db_conn)
    jump_table = SvJumpTable(db_conn) # initialize jump table
    param = ParameterSetManager()
    param.set_selected("SV-PacBio")
    pack, _, _, ret_query_genome = load_genomes(query_genome, reference_genome, param)


    seeds_list = call_table.calls_to_seeds_by_id_auto(pack, run_ids, True, 0)

    reconstructed_query_genome = call_table.reconstruct_sequenced_genome_from_seeds(seeds_list, pack)
    reconstructed_query_genome.store(reconstructed_query_genome_path + "/ma/genome")

    #seeds_n_rects_reconstr = compute_seeds(reconstructed_query_genome_path, query_genome, db_name, 1)


    seeds_list_display = [(reconstructed_query_genome.start_of_sequence(name), seeds, [], []) for name, seeds, _ in seeds_list]

    seeds_n_rects = compute_seeds(query_genome, reference_genome, db_name, 1)

    out = []
    out.append(render_seeds(seeds_list_display, reconstructed_query_genome_path, reference_genome,
                              "reconstructed on reference", "Reconstructed Genome", "Reference Genome"))
    out.append(render_seeds(seeds_n_rects, query_genome, reference_genome, "assembly on reference",
                            "Sequenced Genome", "Reference Genome"))
    #out.append(render_seeds(seeds_n_rects_reconstr, reconstructed_query_genome_path, query_genome,
    #                          "reconstructed on assembly"))

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
                        y += l
                    else:
                        x += l
                xs.append(x + x_start)
                ys.append(y + y_start)
            xs.append(float("NaN"))
            ys.append(float("NaN"))
            iden = 100 * matches / l_total
            print(name, nw_alignment.get_score(), matches, mismatches, indels, indelops, iden, sep="\t")
            break
        plot = figure(title="alignments reconstructed on assembly", plot_width=1000, plot_height=1000)
        plot.xaxis.axis_label = "Sequenced Genome"
        plot.yaxis.axis_label = "Reconstructed Genome"
        for x in cx:
            plot.line(y=[x, x], x=[0, cy[-1]], color="black")
        for y in cy:
            plot.line(y=[0, cx[-1]], x=[y, y], color="black")
        plot.line(x=xs, y=ys, line_width=4)
        out.append(plot)

    show(row(out))

"""
name    score       matches missmatches     indels  indel ops       % identity
chr1    434082      219848  0       1247    920     99.58417509942655
chr2    299750      574595  85160   142053  56703   78.5716141506712
chr3    608894      308063  0       1206    1205    99.61037029372584
chr4    -1208226    817711  379702  384677  141983  58.75973147803492
chr5    1135982     577393  0       3297    3118    99.46288700817038
chr6    106         153     4       282800  13      0.05407182009987383
chr7    -517024     177887  93499   846960  33576   16.08429546270368
chr8    -860        377     48      546241  69      0.06896400335124812
chr9    -801522     147732  217448  754041  14984   13.311875549660472
chr10   -1161596    219792  269995  320413  71505   29.26139578476781
chr11   -653602     219154  139554  464473  40936   27.19636246759516
chr12   -1305272    347611  349548  406525  52381   32.53526492096221
chr13   -43778      17169   3191    850231  3139    1.9739837473153636
chr14   -743488     511724  176595  304046  116035  58.67351560215787
chr15   -1052266    673227  247627  424202  152199  56.02554169582813
chr16   -844228     567840  190310  377936  127638  55.16164097721901
"""
