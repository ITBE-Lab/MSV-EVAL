from yeast.display_mems import *


#query_genome = "knowlesiStrain"
#query_genome = "UFRJ50816-chrVII-section"
query_genome = genome_dir + "UFRJ50816"
reconstructed_query_genome_path = "/MAdata/genome/reconstructed/yeast/UFRJ50816"
#query_genome = "YPS138-chrVII-section"
#reference_genome = "YPS138-chrVII-section"
reference_genome = genome_dir + "YPS138"
#reference_genome = "vivax"

def nw_comparison(reconstructed_query_genome, ret_query_genome,
                  title="alignments reconstructed on assembly"):
    print("name", "score", "% of max score", "matches", "missmatches", "indels", "indel ops", "% identity",
            sep="\t")
    xs = []
    cx = [0]
    ys = []
    cy = [0]
    # write cigar memory dumps to ssd
    MS.util.ksw_file_prefix = "/MAdata/tmp/.CIGARMemoryManager"
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
        param.by_name("Minimal Bandwidth in Gaps").set(10000)
        aligner = NWAlignment(param)
        nw_alignment = aligner.execute(reconstr, assembly)
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

        percent_of_max_score = (100.0 * nw_alignment.get_score()) / (param.by_name("Match Score").get() * l_total)

        print(name, nw_alignment.get_score(), percent_of_max_score, matches, mismatches, indels, indelops, iden,
                sep="\t")
    plot = figure(title=title, plot_width=1000, plot_height=1000)
    decorate_plot(plot, reconstructed_query_genome_path, query_genome)
    plot.xaxis.axis_label = "Sequenced Genome"
    plot.yaxis.axis_label = "Reconstructed Genome"
    plot.line(x=xs, y=ys, line_width=4)
    return plot

if __name__ == "__main__":
    db_name = "UFRJ50816"
    run_ids = [1, 2, 3, 4]
    #run_ids = [5, 6]

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


    print("reconstructing...")
    seeds_list = call_table.calls_to_seeds_by_id(pack, run_ids, True, True)

    reconstructed_query_genome = call_table.reconstruct_sequenced_genome(seeds_list, pack)
    print("done")
    reconstructed_query_genome.store(reconstructed_query_genome_path + "/ma/genome")

    out = []
    print("contig_name", "nt in seeds", "nt in insertions", "nt in ends", "sequenced len", "reconstr len", sep="\t")
    buckets = {}
    nt_insertion_total = 0
    nt_seeds_total = 0
    nt_ends_total = 0
    req_len_total = 0
    num_insertions = 0
    num_insertions_size_one = 0
    num_insertions_larger_eq_hundred = 0
    for idx, ((contig_name, seeds, insertions), reconstr_contig, (y_start, assembly)) in enumerate(zip(seeds_list,
                                                  reconstructed_query_genome.contigNucSeqs(),
                                                  ret_query_genome)):
        nt_in_seeds = 0
        nt_in_ins = 0
        if len(insertions[0]) > 0:
            num_insertions += 1
        if len(insertions[0]) == 1:
            num_insertions_size_one += 1
        if len(insertions[0]) >= 100:
            num_insertions_larger_eq_hundred += 1
        if len(insertions[-1]) > 0:
            num_insertions += 1
        if len(insertions[-1]) == 1:
            num_insertions_size_one += 1
        if len(insertions[-1]) >= 100:
            num_insertions_larger_eq_hundred += 1
        nt_in_ends = len(insertions[0]) + len(insertions[-1])
        nt_ends_total += nt_in_ends
        ins_len = []
        for seed in seeds:
            nt_in_seeds += seed.size
            nt_seeds_total += seed.size
        for nuc_seq in insertions[1:-1]:
            nt_in_ins += len(nuc_seq)
            ins_len.append(len(nuc_seq))
            nt_insertion_total += len(nuc_seq)
            if len(nuc_seq) > 0:
                num_insertions += 1
            if len(nuc_seq) == 1:
                num_insertions_size_one += 1
            if len(nuc_seq) >= 100:
                num_insertions_larger_eq_hundred += 1
        req_len_total += len(reconstr_contig)
        print(contig_name, nt_in_seeds, nt_in_ins, nt_in_ends, len(assembly), len(reconstr_contig), sep="\t")
        for l in ins_len:
            if not l in buckets:
                buckets[l] = 0
            buckets[l] += 1
    print("nt in insertion total:", nt_insertion_total)
    print("nt in seeds total:", nt_seeds_total)
    print("nt in ends total:", nt_ends_total)
    print("reconstr len total:", req_len_total)
    if nt_insertion_total + nt_seeds_total + nt_ends_total != req_len_total:
        print("WARNING: lengths do not match up:", nt_insertion_total + nt_seeds_total)
    print("there are", num_insertions, "insertions in total,", num_insertions_size_one, "are of size 1 nt and",
          num_insertions_larger_eq_hundred, "are of size >= 100 nt. Full distribution shown in plot.")
    if False:
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
    if False:
        out.append(render_seeds(seeds_n_rects, query_genome, reference_genome, "assembly on reference",
                                "Sequenced Genome", "Reference Genome"))

    seeds_list_display = [(reconstructed_query_genome.start_of_sequence(name), seeds, [], []) for name, seeds, _ in seeds_list]
    if False:
        out.append(render_seeds(seeds_list_display, reconstructed_query_genome_path, reference_genome,
                                "reconstructed on reference", "Reconstructed Genome", "Reference Genome"))
    if True:
        out.append(render_seeds(seeds_list_display, reconstructed_query_genome_path, reference_genome,
                                "reconstructed on reference (x & y squeezed)", "Reconstructed Genome",
                                "Reference Genome", True))

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
    if False:
        print("NW comparison for reconstruction & sequenced genome")
        plot = nw_comparison(reconstructed_query_genome, ret_query_genome)
        out.append(plot)

    if False:
        print("NW comparison for reference & sequenced genome")
        plot = nw_comparison(pack, ret_query_genome, title="alignments reference on assembly")
        out.append(plot)

    show(row(out))
