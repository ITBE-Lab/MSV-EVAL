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
    run_ids = [1, 2, 3, 4]
    #run_ids = [1, 2]

    db_conn = DbConn({"SCHEMA": {"NAME": db_name}})
    call_table = SvCallTable(db_conn)

    #call_table.copy_path(3, 2, 110)
    #call_table.copy_path(5, 1, 0)

    jump_table = SvJumpTable(db_conn) # initialize jump table
    param = ParameterSetManager()
    param.set_selected("SV-PacBio")
    pack, _, _, ret_query_genome = load_genomes(query_genome, reference_genome, param)


    seeds_list = call_table.calls_to_seeds_by_id(pack, run_ids, True, 0)

    reconstructed_query_genome = call_table.reconstruct_sequenced_genome(seeds_list, pack)
    reconstructed_query_genome.store(reconstructed_query_genome_path + "/ma/genome")

    out = []
    print("contig_name", "nt in seeds", "nt in insertions", "nt in ends", "sequenced len", "reconstr len", sep="\t")
    buckets = {}
    for idx, ((contig_name, seeds, insertions), reconstr_contig, (y_start, assembly)) in enumerate(zip(seeds_list,
                                                  reconstructed_query_genome.contigNucSeqs(),
                                                  ret_query_genome)):
        nt_in_seeds = 0
        nt_in_ins = 0
        nt_in_ends = len(insertions[0]) + len(insertions[-1])
        ins_len = []
        for seed in seeds:
            nt_in_seeds += seed.size
        for nuc_seq in insertions[1:-1]:
            nt_in_ins += len(nuc_seq)
            ins_len.append(len(nuc_seq))
        print(contig_name, nt_in_seeds, nt_in_ins, nt_in_ends, len(assembly), len(reconstr_contig), sep="\t")
        for l in ins_len:
            if not l in buckets:
                buckets[l] = 0
            buckets[l] += 1
    indel_distrib = figure(title="Insertion distrib", y_axis_type="log", plot_width=1000, plot_height=1000)
    indel_distrib.vbar(x=[key for key, _ in buckets.items()],
                        width=4/5,
                        top=[val for _, val in buckets.items()],
                        bottom=0.1,
                        color="blue")
    indel_distrib.xaxis.axis_label = "Amount"
    indel_distrib.yaxis.axis_label = "Insertion length"
    out.append(indel_distrib)


    seeds_n_rects = compute_seeds(query_genome, reference_genome, db_name, 1)
    out.append(render_seeds(seeds_n_rects, query_genome, reference_genome, "assembly on reference",
                            "Sequenced Genome", "Reference Genome"))

    seeds_list_display = [(reconstructed_query_genome.start_of_sequence(name), seeds, [], []) for name, seeds, _ in seeds_list]
    out.append(render_seeds(seeds_list_display, reconstructed_query_genome_path, reference_genome,
                              "reconstructed on reference", "Reconstructed Genome", "Reference Genome"))

    if False:
        seeds_n_rects_reconstr = compute_seeds(reconstructed_query_genome_path, query_genome, db_name, 1)
        out.append(render_seeds(seeds_n_rects_reconstr, reconstructed_query_genome_path, query_genome,
                                "reconstructed on assembly", "Reconstructed Genome", "Sequenced Genome"))

    if False: # exact match comparison
        print("name", "perfect match", sep="\t")
        for name, reconstr, (y_start, assembly) in zip(
                                                reconstructed_query_genome.contigNames(),
                                                reconstructed_query_genome.contigNucSeqs(),
                                                ret_query_genome):
            print(name, reconstr.equals(assembly), sep="\t")
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
            printed_error = False
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
                    if not printed_error:
                        printed_error = True
                        print("Mismatch:", x, y, reconstr[y-l], "!=", assembly[x-l])
                if op == MatchType.insertion or op == MatchType.deletion:
                    indels += l
                    indelops += 1
                    if op == MatchType.insertion:
                        y += l
                    else:
                        x += l
                    if not printed_error:
                        printed_error = True
                        print("Indel:", x, y)
                xs.append(x + x_start)
                ys.append(y + y_start)
            xs.append(float("NaN"))
            ys.append(float("NaN"))
            iden = 100 * matches / l_total
            print(name, nw_alignment.get_score(), matches, mismatches, indels, indelops, iden, sep="\t")
        plot = figure(title="alignments reconstructed on assembly", plot_width=1000, plot_height=1000)
        decorate_plot(plot, reconstructed_query_genome_path, query_genome)
        plot.xaxis.axis_label = "Sequenced Genome"
        plot.yaxis.axis_label = "Reconstructed Genome"
        plot.line(x=xs, y=ys, line_width=4)
        out.append(plot)

    show(row(out))

"""
contig_name     nt in seeds     nt in insertions        nt in ends      sequenced len   reconstr len
chr1            205485             14038                654             220177          220177
chr2            696721             30811                3769            731301          731301
chr3            304050             3621                 393             308064          308064
chr4            1338350            45047                8221            1391618         1391618
chr5            551437             25561                574             577572          577572
chr6            277427             4826                 704             282957          282957
chr7            1053727            51673                567             1105967         1105967
chr8            519975             26251                436             546662          546662
chr9            328747             45315                563             374625          374625
chr10           689704             40243                21186           751133          751133
chr11           789360             15916                545             805821          805821
chr12           1001039            59147                8227            1068413         1068413
chr13           811600             49844                8320            869764          869764
chr14           852128             19503                524             872155          872155
chr15           1163943            28849                8851            1201643         1201643
chr16           1009782            19463                166             1029411         1029411
"""
