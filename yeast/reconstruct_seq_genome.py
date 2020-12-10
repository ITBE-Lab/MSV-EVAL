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
    #run_ids = [1, 2, 3, 4]
    run_ids = [5, 6]

    db_conn = DbConn({"SCHEMA": {"NAME": db_name}})
    call_table = SvCallTable(db_conn)

    # copy order from 1 to 5 (pacBio)
    call_table.copy_path(1, 5, 100)
    # copy order from 3 to 6 (Illumina)
    call_table.copy_path(3, 6, 0)

    jump_table = SvJumpTable(db_conn) # initialize jump table
    param = ParameterSetManager()
    param.set_selected("SV-PacBio")
    pack, _, _, ret_query_genome = load_genomes(query_genome, reference_genome, param)


    seeds_list = call_table.calls_to_seeds_by_id(pack, run_ids, True, 0)

    print("reconstructing...")
    reconstructed_query_genome = call_table.reconstruct_sequenced_genome(seeds_list, pack)
    print("done")
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
        print("name", "score", "% of max score", "matches", "missmatches", "indels", "indel ops", "% identity",
                sep="\t")
        xs = []
        cx = [0]
        ys = []
        cy = [0]
        for y_start, name, reconstr, (x_start, assembly) in zip(
                                                reconstructed_query_genome.contigStarts(),
                                                reconstructed_query_genome.contigNames(),
                                                reconstructed_query_genome.contigNucSeqs(),
                                                ret_query_genome):
            x = 0
            y = 0
            xs.append(x_start)
            ys.append(y_start)
            cx.append(x_start + len(assembly))
            cy.append(y_start + len(reconstr))
            matches = 0
            mismatches = 0
            indels = 0
            indelops = 0
            l_total = min(len(reconstr), len(assembly))
            nw_alignment = runKsw(reconstr, assembly, 10000)
            printed_error = True # set to false in order to print the postion of the first error
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

            match_score = param.by_name("Match Score").get()
            gap_penalty = 0
            if len(reconstr) != len(assembly):
                gap_open_penalty = param.by_name("Gap penalty").get()
                gap_open_penalty_2 = param.by_name("Second Gap penalty").get()
                gap_extend_penalty = param.by_name("Extend Penalty").get()
                gap_extend_penalty_2 = param.by_name("Second Extend Penalty").get()
                gapLen = abs(len(reconstr) - len(assembly))
                gap_penalty = min(gap_open_penalty + gap_extend_penalty * gapLen,
                                  gap_open_penalty_2 + gap_extend_penalty_2 * gapLen)
            percent_of_max_score = (100.0 * nw_alignment.get_score()) / (match_score * l_total - gap_penalty)

            print(name, nw_alignment.get_score(), percent_of_max_score, matches, mismatches, indels, indelops, iden,
                  sep="\t")
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


"""
chr1    150029  11825   1       220177  161855
chr2    348831  21451   1       731301  370283
chr3    171755  2434    0       308064  174189
chr4    704895  31696   2       1391618 736593
chr5    276762  17546   2       577572  294310
chr6    179867  1809    1       282957  181677
chr7    794300  50328   0       1105967 844628
chr8    312796  22533   1       546662  335330
chr9    214000  9502    2       374625  223504
chr10   439853  35692   5       751133  475550
chr11   426265  11321   1       805821  437587
chr12   600929  27645   1       1068413 628575
chr13   1153526 46186   2       869764  1199714
chr14   460195  17037   1       872155  477233
chr15   597630  23566   2       1201643 621198
chr16   532746  17050   2       1029411 549798
name    score   % of max score  matches missmatches     indels  indel ops       % identity
chr1    433574  99.86870715707231       218253  324     9386    86      99.12615759139238
chr2    1430666 99.25813069471242       718865  486     16827   202     99.25948734376468
chr3    608188  99.64169567888592       305622  341     9931    76      99.20730757245248
chr4    2728520 100.45797603237901      1367793 408     24457   145     99.89424798118081
chr5    1137284 99.58015163574498       572040  613     5491    166     99.7932748920581
chr6    527774  98.0852182866021        265665  292     24729   87      97.06926916247086
chr7    2193042 99.992157628029         1098983 292     7152    125     99.93161989024628
chr8    1074774 98.61669037023444       540218  570     15198   131     98.82120944934896
chr9    678316  102.32584801885056      342029  644     35128   143     98.89547172320869
chr10   1426034 95.94949382129597       740585  4858    27388   3276    98.59572139687646
chr11   1596412 99.47738003326273       801502  525     5315    178     99.74537924305704
chr12   2096408 100.25513217376916      1051670 586     17071   160     99.85757285148647
chr13   1687042 100.71904563468135      846861  719     22866   185     99.83483876443835
chr14   1732516 99.65264256810973       869303  588     10259   159     99.67299390589976
chr15   2357332 99.74181635704657       1184944 1288    17546   417     99.71195766964246
chr16   2030518 100.11971855292418      1018132 516     11290   128     99.89766232491966

"""