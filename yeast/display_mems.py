from bokeh.plotting import figure, show, reset_output, ColumnDataSource
from bokeh.layouts import column, row
from bokeh.models import FuncTickFormatter
from bokeh.models.tools import HoverTool
from bokeh.io import output_file
from MA import *
from MSV import *
from sv_util.os_aligners import *
from bisect import bisect_right

global_prefix = "/MAdata/"

def load_genomes(query_genome, reference_genome, params):
    pack = Pack()
    pack.load(reference_genome + "/ma/genome")
    fm_index = FMIndex()
    #fm_index.load(reference_genome + "/ma/genome") # not used at the moment
    mm_index = MinimizerIndex(params, libMA.util.StringVector([str(pack.extract_forward_strand())]), libMA.util.StringVector(["chrVII"]))
    file_reader = FileReader(params)
    ret_query_genome = []
    if False: # load query genome from fasta
        f_stream = FileStreamFromPath(query_genome + "/fasta/genome.fna")
        idx = 0
        y_start = 0
        while not f_stream.eof():
            ret_query_genome.append((y_start, file_reader.execute(f_stream)))
            y_start += len(ret_query_genome[-1][1])
            ret_query_genome[-1][1].name = "sequence" + str(idx)
            idx += 1
    else: # load query genome from pack
        query_pack = Pack()
        query_pack.load(query_genome + "/ma/genome")
        ret_query_genome = list(zip(query_pack.contigStarts(), query_pack.contigNucSeqs()))


    #print("seq", pack.extract_from_to(185500, 190900))
    #exit()

    return pack, fm_index, mm_index, ret_query_genome

def compute_seeds(query_genome, reference_genome, db_name, seq_id, ambiguity=2):
    param = ParameterSetManager()
    param.set_selected("SV-PacBio")

    pack, fm_index, mm_index, query_genomes = load_genomes(query_genome, reference_genome, param)
    mm_index.set_max_occ(ambiguity)

    seeding_module = MMFilteredSeeding(param)
    seed_lumper = SeedLumping(param)
    soc_module = StripOfConsiderationSeeds(param)
    soc_filter = GetAllFeasibleSoCsAsSet(param)
    reseeding_1 = RecursiveReseedingSoCs(param, pack)
    db_conn = DbConn({"SCHEMA": {"NAME": db_name}})
    mm_counter = HashFilterTable(db_conn).get_counter(seq_id)
    jumps_from_seeds = SvJumpsFromSeeds(param, pack)

    ret = []
    for y_start, query_genome in query_genomes:
        seeds = seeding_module.execute(mm_index, query_genome, pack, mm_counter)
        mems = seed_lumper.execute(seeds, query_genome, pack)
        #ret.append((y_start, mems, []))
        #continue
        socs = soc_module.execute(mems, query_genome, pack)
        soc_filtered_seeds = soc_filter.execute(socs)
        filtered_seeds_2, helper_ret_1 = reseeding_1.cpp_module.execute_helper(soc_filtered_seeds, pack, query_genome)

        filtered_seeds = []
        #while not socs.empty():
        #    for seed in socs.pop():
        #        filtered_seeds.append( (seed, "initial SoC") )
        #for seed in helper_ret_1.seed_removed:
        #    filtered_seeds.append( (seed, "reseeding Soc / overlapping filter / enclosed SoC") )
        #for seed, parlindrome, overlapping in zip(helper_ret_1.seeds, helper_ret_1.parlindrome,
        #                                          helper_ret_1.overlapping):
        #    if parlindrome:
        #        filtered_seeds.append((seed, "palrindrome"))
        #    if overlapping:
        #        filtered_seeds.append((seed, "overlapping"))


        #r = []
        #for seeds in soc_filtered_seeds.content:
        #    for seed in seeds:
        #        r.append(seed)

        ret.append((y_start, filtered_seeds_2, helper_ret_1.rectangles, filtered_seeds))

        #ret.append((y_start, r, helper_ret_1.rectangles, filtered_seeds))
    return ret


def run_aligner(query_genome_str, reference_genome):
    pack, fm_index, mm_index, query_genomes = load_genomes(query_genome_str, reference_genome)

    read_set = {"technology":"pb", "name":"test", "fasta_file":"sequenced.fasta"}
    json_dict = {"reference_path":reference_genome}
    sam_file_path = "mm2.sam"
    ret = []
    for y_start, query_genome in query_genomes:
        with open("sequenced.fasta", "w") as out_file:
            out_file.write(">sequence0\n")
            out_file.write(str(query_genome))
            out_file.write("\n")

        mm2(read_set, sam_file_path, json_dict)
        #ngmlr(read_set, sam_file_path, json_dict)
        read = ReadByName()
        query_genome.name = "sequence0"
        read.append(query_genome)

        f_stream = FileStreamFromPath("mm2.sam")
        file_reader = SamFileReader(ParameterSetManager())
        while not f_stream.eof():
            alignment = file_reader.execute(f_stream, pack, read)
            ret.append((y_start, alignment))

    return ret

def run_aligner_seeds(query_genome_str, reference_genome):
    ret = []
    for y_start, alignment in run_aligner(query_genome_str, reference_genome):
        seeds_list = []
        seeds = alignment.to_seeds(pack)
        for seed in seeds:
            seeds_list.append(seed)

        ret.append((y_start, seeds_list))

    return ret

def seeds_to_jumps(seeds_n_rects, query_genome, reference_genome):
    param = ParameterSetManager()
    # make sure we record the entire genome
    param.by_name("Maximal Dummy Distance").set(10000000)
    param.by_name("Minimal Dummy Distance").set(0)

    pack, fm_index, mm_index, query_genomes = load_genomes(query_genome, reference_genome, param)
    jumps_from_seeds = SvJumpsFromExtractedSeeds(param, pack)

    ret = []
    for (y_start, seeds, _, _), (_, query_genome) in zip(seeds_n_rects, query_genomes):
        jumps = jumps_from_seeds.execute(libMA.containers.Seeds(seeds), pack, query_genome)
        for jump in jumps:
            jump.query_from += y_start
            jump.query_to += y_start
        ret.append(jumps)
    return ret

def decorate_plot(plot, query_genome, reference_genome, diagonal=False, x_only=False, max_y=None):
    if not x_only:
        pack_1 = Pack()
        pack_1.load(query_genome + "/ma/genome")
    pack_2 = Pack()
    pack_2.load(reference_genome + "/ma/genome")

    xs = []
    ys = []
    for idx in [*pack_2.contigStarts(), pack_2.unpacked_size_single_strand]:
        xs.append(idx)
        ys.append(0)
        xs.append(idx)
        if not x_only:
            ys.append(pack_1.unpacked_size_single_strand)
        else:
            ys.append(max_y)
        xs.append(float("NaN"))
        ys.append(float("NaN"))
    plot.line(x=xs, y=ys, color="black", line_width=1)
    if not x_only:
        ys = []
        xs = []
        for idx in [*pack_1.contigStarts(), pack_1.unpacked_size_single_strand]:
            xs.append(0)
            ys.append(idx)
            xs.append(pack_2.unpacked_size_single_strand)
            ys.append(idx)
            xs.append(float("NaN"))
            ys.append(float("NaN"))
        plot.line(x=xs, y=ys, color="black", line_width=1)

    if diagonal:
        plot.line(x=[0, pack_2.unpacked_size_single_strand], y=[0, pack_1.unpacked_size_single_strand],
                  color="black", line_width=1)

    ticker_code = """
            if(tick < 0 || tick >= genome_end)
                return "n/a";
            idx = 0;
            while(contig_starts[idx + 1] < tick)
                idx += 1;
            return contig_names[idx] + ": " + (tick - contig_starts[idx]);
        """
    plot.xaxis[0].formatter = FuncTickFormatter(
                    args={"contig_starts": [*pack_2.contigStarts(), pack_2.unpacked_size_single_strand],
                            "genome_end":pack_2.unpacked_size_single_strand,
                            "contig_names": [*pack_2.contigNames()]},
                    code=ticker_code)
                    
    plot.xaxis.major_label_orientation = math.pi/4
    if not x_only:
        plot.yaxis[0].formatter = FuncTickFormatter(
                        args={"contig_starts": [*pack_1.contigStarts(), pack_1.unpacked_size_single_strand],
                                "genome_end":pack_1.unpacked_size_single_strand,
                                "contig_names": [*pack_1.contigNames()]},
                        code=ticker_code)

def render_jumps(jumps, query_genome, reference_genome):
    xs = []
    ys = []
    cs = []
    ms = []
    ids = []
    qlens = []
    fs = []
    ts = []
    pack_2 = Pack()
    pack_2.load(reference_genome + "/ma/genome")
    contig_filter = JumpsFilterContigBorder(ParameterSetManager())
    colors = {True:{True:"blue", False:"green"}, False:{True:"purple", False:"orange"}}
    for idx, jumps_ in enumerate(jumps):
        for jump in jumps_:
            #if idx != 8281:
            #    idx+=1
            #    continue
            f = jump.from_pos
            t = jump.to_pos
            if not jump.from_known():
                f = t
            if not jump.to_known():
                t = f
            #if abs(f - t) < 200:
            #    continue
            xs.append(f)
            ys.append(t)
            ms.append(jump.was_mirrored)
            if contig_filter.cpp_module.by_contig_border(jump, pack_2):
                cs.append("red")
            else:
                cs.append(colors[jump.from_forward][jump.to_forward])
            ids.append(jump.id)
            qlens.append( jump.query_to - jump.query_from )
            fs.append("fow" if jump.from_forward else "rev")
            ts.append("fow" if jump.to_forward else "rev")
    plot = figure(title="entries", plot_width=1000, plot_height=1000)
    decorate_plot(plot, reference_genome, reference_genome, True)
    plot.x(x="xs", y="ys", color="cs", line_width=4, source=ColumnDataSource(data={
                    "xs":xs, "ys":ys, "cs":cs, "ms":ms, "ids": ids, "qlens": qlens, "fs":fs, "ts":ts,
                }))
    plot.add_tools(HoverTool(tooltips=[("x, y", "@xs, @ys"),("id", "@ids"),("color, mirrored", "@cs, @ms"),
                                       ("query dist", "@qlens"), ("strandinfo", "from @fs to @ts")]))
    plot.xaxis.axis_label = "Reference Genome"
    plot.yaxis.axis_label = "Reference Genome"
    return plot

def jumps_to_calls_to_db(jumps, db_name, query_genome_str, reference_genome, min_cov=10):
    parameter_set_manager = ParameterSetManager()
    parameter_set_manager.set_selected("SV-PacBio")
    min_size = parameter_set_manager.by_name("Min Size Edge").get()
    db_conn = DbConn({"SCHEMA": {"NAME": db_name}})

    JumpRunTable(db_conn)
    SvCallerRunTable(db_conn)
    sv_call_table = SvCallTable(db_conn)
    #call_per_contig_table = FirstCallPerContigTable(db_conn)
    caller_run_id = GetCallInserter(parameter_set_manager, db_conn, "Ground Truth - Large",
                                   "SV's form genome assembly that are too long for Illumina reads", -1).cpp_module.id
    contig_border_caller_run_id = GetCallInserter(parameter_set_manager, db_conn, "Ground Truth - Contig Border",
                                   "SV's form genome assembly that are not covered by any alignments", -1).cpp_module.id
    small_caller_run_id = GetCallInserter(parameter_set_manager, db_conn, "Ground Truth - Small",
                                   "SV's form genome assembly that are too small for pacBio reads", -1).cpp_module.id
    contig_start_caller_run_id = GetCallInserter(parameter_set_manager, db_conn, "Ground Truth - Contig Start",
                                   "SV's form genome assembly that are not covered by any alignments", -1).cpp_module.id

    pack, fm_index, mm_index, query_genomes = load_genomes(query_genome_str, reference_genome, parameter_set_manager)


    contig_filter = JumpsFilterContigBorder(parameter_set_manager)

    cnt_small = 0
    cnt_large = 0
    cnt_contig_border = 0

    for contig_idx, (jumps_, (y_start, query_genome)) in enumerate(zip(jumps, query_genomes)):
        query_genome.check()
        is_first = True
        idx = 0
        for jump in jumps_:
            jump.id = idx
            assert(jump.num_supp_nt() > 0)
            f = jump.from_pos
            t = jump.to_pos
            if not jump.from_known():
                if jump.from_forward:
                    f = pack.start_of_sequence_id(pack.seq_id_for_pos(t))
                else:
                    f = pack.end_of_sequence_id(pack.seq_id_for_pos(t))
            if not jump.to_known():
                if jump.from_forward:
                    t = pack.end_of_sequence_id(pack.seq_id_for_pos(f))
                else:
                    t = pack.start_of_sequence_id(pack.seq_id_for_pos(f))

            call = SvCall(f, t, 0, 0, jump.from_forward, jump.to_forward, 1, jump.num_supp_nt())
            if jump.query_from < jump.query_to:
                call.inserted_sequence = NucSeq(query_genome, jump.query_from - y_start, jump.query_to - y_start)
                call.inserted_sequence.check()
            call.order_id = idx
            idx += 1
            call.ctg_order_id = contig_idx
            call.mirrored = jump.was_mirrored
            from_forward = jump.from_forward
            if jump.was_mirrored:
                from_forward = not jump.to_forward
            if not jump.from_known() or not jump.to_known():
                sv_call_table.insert_call(contig_start_caller_run_id, call)
            elif contig_filter.cpp_module.by_contig_border(jump, pack):
                cnt_contig_border += 1
                sv_call_table.insert_call(contig_border_caller_run_id, call)
            elif max(abs(f - t), abs(jump.query_from - jump.query_to)) < min_size:
                cnt_small += 1
                sv_call_table.insert_call(small_caller_run_id, call)
            else:
                cnt_large += 1
                sv_call_table.insert_call(caller_run_id, call)
            #if is_first:
            #    is_first = False
            #    call_per_contig_table.insert(call.id, query_genome.name, from_forward)

    print("Inserted into DB. There are", cnt_small, "small and", cnt_large, "large entries.", cnt_contig_border,
          "entries are too close to contig borders and are filtered out.")

    return caller_run_id, None

def render_seeds(seeds_1, query_genome, reference_genome, title="seeds", y_axis="Sequenced Genome",
                 x_axis="Reference Genome"):
    plot = figure(title=title, plot_width=1000, plot_height=1000)

    decorate_plot(plot, query_genome, reference_genome)

    if True:
        xs = []
        xe = []
        ys = []
        ye = []
        for y_start, _, rects, _ in seeds_1:
            for rect in rects:
                xs.append(rect.x_axis.start)
                ys.append(rect.y_axis.start + y_start)
                xe.append(rect.x_axis.start + rect.x_axis.size)
                ye.append(rect.y_axis.start + rect.y_axis.size + y_start)
        plot.quad(left="xs", bottom="ys", right="xe", top="ye", fill_color="black",
                        fill_alpha=0.2, line_width=0,
                        source=ColumnDataSource({"xs":xs, "xe":xe, "ys":ys, "ye":ye}))

    xs = {True:{True:[], False:[]}, False:{True:[], False:[]}}
    ys = {True:{True:[], False:[]}, False:{True:[], False:[]}}
    rs = {True:{True:[], False:[]}, False:{True:[], False:[]}}
    cs = {True:{True:"purple", False:"green"}, False:{True:"blue", False:"orange"}}
    lw = {True:3, False: 4}
    for y_start, seeds, _, filtered_seeds in seeds_1:
        for seedlist, filtered in [(seeds, False), (filtered_seeds, True)]:
            for seed_ in seedlist:
                if filtered:
                    # uaaaagh @todo better construct needed
                    seed, reason = seed_
                else:
                    seed = seed_
                    reason = None
                #if seed.size < 20:
                #    continue
                xs[filtered][seed.on_forward_strand].append([seed.start_ref,
                                            seed.start_ref + seed.size * (1 if seed.on_forward_strand else -1)])
                ys[filtered][seed.on_forward_strand].append([seed.start + y_start, seed.start + seed.size + y_start])
                rs[filtered][seed.on_forward_strand].append(reason)
    for filtered in [True, False]:
        for forw in [True, False]:
            plot.multi_line(xs="xs", ys="ys", color=cs[filtered][forw], line_width=lw[filtered], line_cap="round",
                            source=ColumnDataSource(data={
                                "xs":xs[filtered][forw], "ys":ys[filtered][forw], "rs":rs[filtered][forw]
                            }))

    plot.xaxis.axis_label = x_axis
    plot.yaxis.axis_label = y_axis
    plot.add_tools(HoverTool(tooltips=[("filtered due to", "@rs")]))
    return plot

genome_dir = global_prefix + "genome/yeasts/"

#query_genome = "knowlesiStrain"
#query_genome = "UFRJ50816-chrVII-section"
query_genome = genome_dir + "UFRJ50816"
#query_genome = genome_dir + "UWOPS919171"
#query_genome = "YPS138-chrVII-section"
#reference_genome = "YPS138-chrVII-section"
reference_genome = genome_dir + "YPS138"
#reference_genome = genome_dir + "SK1"
#reference_genome = "vivax"

db_name = "UFRJ50816"
seq_id = 1

if __name__ == "__main__":
    out = []

    seeds_n_rects = compute_seeds(query_genome, reference_genome, db_name, seq_id)
    jumps = seeds_to_jumps(seeds_n_rects, query_genome, reference_genome)

    jumps_to_calls_to_db(jumps, db_name, query_genome, reference_genome)
    #exit()

    out.append(render_seeds(seeds_n_rects, query_genome, reference_genome))
    out.append(render_jumps(jumps, query_genome, reference_genome))
    show(row(out))


# ~/workspace/samtools/samtools view -h minimap2.sorted.bam "chrIX:251682-259682" > minimap2.filtered.bam


# ~/workspace/samtools/samtools view minimap2.filtered.bam | awk '{print length($10), $1}' | sort -n
