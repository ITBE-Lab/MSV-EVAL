from sv_util.settings import *
from MSV import *
from yeast.load_genomes import *
from bokeh.plotting import figure, show, save

genome_dir = main_data_folder + "/genome/yeasts/"
#query_genome = "knowlesiStrain"
#query_genome = "UFRJ50816-chrVII-section"
query_genome = genome_dir + "UFRJ50816"
#query_genome = "YPS138-chrVII-section"
#reference_genome = "YPS138-chrVII-section"
reference_genome = genome_dir + "YPS138"
#reference_genome = "vivax"

def nw_comparison(reconstructed_query_genome, ret_query_genome,
                  title="alignments reconstructed on assembly", y_axis_label="Reconstructed Genome"):
    print("name", "score", "% of max score", "matches", "missmatches", "indels", "indel ops", "% identity",
            sep="\t")
    xs = []
    cx = [0]
    ys = []
    cy = [0]
    # write cigar memory dumps to ssd
    MS.util.ksw_file_prefix = tmp_ksw_file_prefix
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
    plot.yaxis.axis_label = y_axis_label
    plot.line(x=xs, y=ys, line_width=4)
    return plot

ground_throuth_lt_ten = 18
ground_throuth_ge_ten = 17
ground_throuth_large = 9
to_analyze = [
    #(5, 3, "Gridss"), # Gridss: small, large ## DEFAULT
    (32, 31, 29, 25, 25, "Gridss"), # Gridss: small, large
    #(1, 2, "MA") # MA: small, large
]

if True:
    for lt_ten_id, ge_ten_id, large_id, blur_small, blur_large, name in to_analyze:
        print("Analyzing:", name)
        db_name = "UFRJ50816"
        #run_ids = [3, 4, 5, 6]
        run_ids = [lt_ten_id, ge_ten_id, large_id]

        db_conn = DbConn({"SCHEMA": {"NAME": db_name}})
        call_table = SvCallTable(db_conn)

        print("copying path information... (may take a while)")
        call_table.copy_path(ground_throuth_lt_ten, lt_ten_id, blur_small)
        call_table.copy_path(ground_throuth_ge_ten, ge_ten_id, blur_small)
        call_table.copy_path(ground_throuth_large, large_id, blur_large)
        print("done")

        jump_table = SvJumpTable(db_conn) # initialize jump table
        param = ParameterSetManager()
        param.set_selected("SV-PacBio")
        pack, _, _, ret_query_genome = load_genomes(query_genome, reference_genome, param)


        print("reconstructing...")
        seeds_list = call_table.calls_to_seeds_by_id(pack, run_ids, True, True)

        reconstructed_query_genome = reconstruct_sequenced_genome(seeds_list, pack)
        print("done")
        reconstructed_query_genome.store(reconstructed_query_genome_path + "/ma/genome")

        out = []
        print("caller", "contig_name", "nt in seeds", "nt in insertions", "nt in ends", "sequenced len",
              "reconstr len", sep="\t")
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
            print(name, contig_name, nt_in_seeds, nt_in_ins, nt_in_ends, len(assembly), len(reconstr_contig), sep="\t")
            for l in ins_len:
                if not l in buckets:
                    buckets[l] = 0
                buckets[l] += 1
        print(name, "nt in insertion total:", nt_insertion_total)
        print(name, "nt in seeds total:", nt_seeds_total)
        print(name, "nt in ends total:", nt_ends_total)
        print(name, "reconstr len total:", req_len_total)
        if nt_insertion_total + nt_seeds_total + nt_ends_total != req_len_total:
            print("WARNING: lengths do not match up:", nt_insertion_total + nt_seeds_total)
        print(name, "there are", num_insertions, "insertions in total,", num_insertions_size_one,
             "are of size 1 nt and", num_insertions_larger_eq_hundred, "are of size >= 100 nt.")
        MS.util.ksw_file_system_min_gb_size = ksw_file_system_min_gb_size
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


        if True:
            seeds_n_rects = compute_seeds(query_genome, reference_genome, db_name, 1)
            out.append(render_seeds(seeds_n_rects, query_genome, reference_genome, "assembly on reference",
                                    "Sequenced Genome", "Reference Genome"))

        if True:
            seeds_list_display = [(reconstructed_query_genome.start_of_sequence(name), seeds, [], []) for name, seeds, _ in seeds_list]
            out.append(render_seeds(seeds_list_display, reconstructed_query_genome_path, reference_genome,
                                    "reconstructed on reference", "Reconstructed Genome", "Reference Genome"))
        if False:
            out.append(render_seeds(seeds_list_display, reconstructed_query_genome_path, reference_genome,
                                    "reconstructed on reference (x & y squeezed)", "Reconstructed Genome",
                                    "Reference Genome", True))

        if True:
            seeds_n_rects_reconstr = compute_seeds(reconstructed_query_genome_path, query_genome, db_name, 1)
            out.append(render_seeds(seeds_n_rects_reconstr, reconstructed_query_genome_path, query_genome,
                                    "reconstructed on sequenced genome", "Reconstructed Genome", "Sequenced Genome"))

        if run_ksw:
            print("NW comparison for reconstruction & sequenced genome")
            plot = nw_comparison(reconstructed_query_genome, ret_query_genome)
            out.append(plot)

        if run_ksw:
            print("NW comparison for reference & sequenced genome")
            plot = nw_comparison(pack, ret_query_genome, title="alignments reference on assembly",
                                y_axis_label="Reference Genome")
            out.append(plot)

        output_file(accuracy_recall_data_dir + "/" + name + ".reconstruction.html")
        if show_plots:
            show(row(out))
        if save_plots:
            save(row(out))
