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
    #call_table.copy_path(1, 5, 100)
    # copy order from 3 to 6 (Illumina)
    #call_table.copy_path(3, 6, 0)

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
            gap_open_penalty = param.by_name("Gap penalty").get()
            gap_open_penalty_2 = param.by_name("Second Gap penalty").get()
            gap_extend_penalty = param.by_name("Extend Penalty").get()
            gap_extend_penalty_2 = param.by_name("Second Extend Penalty").get()
            gap_penalty = 0
            if len(reconstr) != len(assembly):
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
name    score   % of max score    matches missmatches     indels  indel ops       % identity
Indel: 2 0
chr1    147752  0       123492  8515    118018  5596    76.29792097865374
using filesystem because: 398825922000 / 96636764160 bytes are required.
that is 371.436 / 90 gigabytes.
^[[BCyclicFileCache hits: 12662144791 misses: 63 = 100%
Indel: 37565 0
chr2    409736  0       322153  47368   362542  1196    87.00183373257751
Indel: 14 0
chr3    318204  0       171881  417     137657  736     98.67500243987853
using filesystem because: 1396207914096 / 96636764160 bytes are required.
that is 1300.32 / 90 gigabytes.
CyclicFileCache hits: 46956358766 misses: 176 = 100%
Indel: 617631 0
chr4    -1700974        0       207962  528592  655103  72      28.232959042510586
using filesystem because: 247865305744 / 96636764160 bytes are required.
that is 230.843 / 90 gigabytes.
CyclicFileCache hits: 7933702195 misses: 36 = 100%
Indel: 32386 0
chr5    515124  0       290112  765     290128  1732    98.57361285719139
Indel: 0 2322
chr6    217898  0       160902  8614    125602  5756    88.56487062203801
using filesystem because: 511773446608 / 96636764160 bytes are required.
that is 476.626 / 90 gigabytes.
CyclicFileCache hits: 25549552388 misses: 86 = 100%
Indel: 107 0
chr7    10420   0       559605  235322  360741  10479   66.25461149760605
using filesystem because: 187292552848 / 96636764160 bytes are required.
that is 174.43 / 90 gigabytes.
CyclicFileCache hits: 7492523883 misses: 25 = 100%
Indel: 3 2
chr8    505084  0       318104  9228    227328  6721    94.86297080487877
Indel: 138 0
chr9    217832  0       200787  13298   169959  10163   89.83597608991337
using filesystem because: 339310054672 / 96636764160 bytes are required.
that is 316.007 / 90 gigabytes.
CyclicFileCache hits: 14052374229 misses: 54 = 100%
Indel: 0 3031
chr10   744078  0       419405  3115    381643  3444    88.19367048680475
using filesystem because: 459145442464 / 96636764160 bytes are required.
that is 427.613 / 90 gigabytes.
CyclicFileCache hits: 15937831046 misses: 75 = 100%
Indel: 25737 0
chr11   175716  0       326492  107045  376334  2143    74.61190574674293
using filesystem because: 748140476784 / 96636764160 bytes are required.
that is 696.76 / 90 gigabytes.
CyclicFileCache hits: 28566848494 misses: 130 = 100%
Indel: 0 1
chr12   -750196 0       297720  330017  441514  752     47.36427633933898
using filesystem because: 684947219568 / 96636764160 bytes are required.
that is 637.907 / 90 gigabytes.
CyclicFileCache hits: 32561515015 misses: 118 = 100%
Indel: 17 0
chr13   -1656092        0       397072  443265  388804  57171   45.652843759916486
using filesystem because: 534292481440 / 96636764160 bytes are required.
that is 497.599 / 90 gigabytes.
CyclicFileCache hits: 18721733612 misses: 89 = 100%
Indel: 15707 0
chr14   -94750  0       307560  168588  397092  1044    64.44650726165206
using filesystem because: 1059930420496 / 96636764160 bytes are required.
that is 987.137 / 90 gigabytes.
CyclicFileCache hits: 34586212516 misses: 131 = 100%
Indel: 526944 0
chr15   -1436462        0       176262  444754  580809  370     28.374527928293393
using filesystem because: 759030533136 / 96636764160 bytes are required.
that is 706.902 / 90 gigabytes.
CyclicFileCache hits: 25813880701 misses: 132 = 100%
Mismatch: 2 2 A != T
chr16   -823002 0       234070  314318  482433  2404    42.5738180204366

"""