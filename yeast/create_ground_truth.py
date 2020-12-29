from bokeh.plotting import figure, show, reset_output, ColumnDataSource
from bokeh.layouts import column, row
from bokeh.models import FuncTickFormatter
from bokeh.models.tools import HoverTool
from bokeh.io import output_file
from MA import *
from MSV import *
from sv_util.os_aligners import *
from sv_util.settings import *
from bisect import bisect_right, bisect_left

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

def decorate_plot_helper(plot, x_axis_start, x_axis_names, y_axis_start, y_axis_names,
                         diagonal=False, x_only=False, max_y=None):
    xs = []
    ys = []
    for idx in x_axis_start:
        xs.append(idx)
        ys.append(0)
        xs.append(idx)
        if not x_only:
            ys.append(y_axis_start[-1])
        else:
            ys.append(max_y)
        xs.append(float("NaN"))
        ys.append(float("NaN"))
    plot.line(x=xs, y=ys, color="black", line_width=1)
    if not x_only:
        ys = []
        xs = []
        for idx in y_axis_start:
            xs.append(0)
            ys.append(idx)
            xs.append(x_axis_start[-1])
            ys.append(idx)
            xs.append(float("NaN"))
            ys.append(float("NaN"))
        plot.line(x=xs, y=ys, color="black", line_width=1)

    if diagonal:
        plot.line(x=[0, x_axis_start[-1]], y=[0, y_axis_start[-1]],
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
                    args={"contig_starts": x_axis_start,
                            "genome_end":x_axis_start[-1],
                            "contig_names": x_axis_names},
                    code=ticker_code)
                    
    plot.xaxis.major_label_orientation = math.pi/4
    if not x_only:
        plot.yaxis[0].formatter = FuncTickFormatter(
                        args={"contig_starts": y_axis_start,
                                "genome_end": y_axis_start[-1],
                                "contig_names": y_axis_names},
                        code=ticker_code)

def decorate_plot(plot, query_genome, reference_genome, diagonal=False, x_only=False, max_y=None):
    if not x_only:
        pack_1 = Pack()
        pack_1.load(query_genome + "/ma/genome")
    pack_2 = Pack()
    pack_2.load(reference_genome + "/ma/genome")
    decorate_plot_helper(plot,
                         [*pack_2.contigStarts(), pack_2.unpacked_size_single_strand],
                         [*pack_2.contigNames()],
                         None if x_only else [*pack_1.contigStarts(), pack_1.unpacked_size_single_strand],
                         None if x_only else [*pack_1.contigNames()],
                         diagonal, x_only, max_y)

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
            if (not jump.from_known() or not jump.to_known()) and not jump.from_forward:
                # flip from and to information for dummy jumps from reverse strand seeds
                # Why do we need to do this?:
                # On reads we do not know the orientation of the reads. So, in order to unify forward and reverse
                # strand reads, we mirror the dummy calls as well. Since they technically do not have a x and y position
                # they cannot get the mirrored flag. Instead we simply mirror all reverse strand dummy calls.
                # Hence these calls are wrong in the context of the reconstruction.
                # since we do not (and cannot since there is no coverage at chrom ends) analyze these calls we can
                # place them as we wish... we choose to simply mirror everything back to the original position
                # since i do not want to do this in the cpp code i do it there in the py code.
                tmp = f
                f = t
                t = tmp

            call = SvCall(f, t, 0, 0, jump.from_forward, jump.to_forward, 1, jump.num_supp_nt())
            if jump.query_from < jump.query_to:
                call.inserted_sequence = NucSeq(query_genome, jump.query_from - y_start, jump.query_to - y_start)
                call.inserted_sequence.check()
            call.order_id = idx
            idx += 1
            call.ctg_order_id = contig_idx
            call.mirrored = jump.was_mirrored
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

class AxisSqueezer():
    def __init__(self, squeeze):
        self.data = set()
        self.squeeze = squeeze
    def add(self, i):
        self.data.add(i)
    def prep(self):
        if self.squeeze:
            self.data = sorted(list(self.data))
    def get(self, i):
        if not self.squeeze:
            return i
        return bisect_left(self.data, i)

def combine_consecutive_diag_seeds(seeds):
    ret = []
    for seed in seeds:
        if len(ret) != 0 and ret[-1].on_forward_strand == seed.on_forward_strand and \
            ((seed.on_forward_strand and seed.start_ref - seed.start == ret[-1].start_ref - ret[-1].start) or \
                (not seed.on_forward_strand and seed.start_ref + seed.start == ret[-1].start_ref + ret[-1].start)):
            ret[-1].size += seed.size
        else:
            ret.append(Seed(seed.start, seed.size, seed.start_ref, seed.on_forward_strand))
    return ret

def render_seeds(seeds_1, query_genome, reference_genome, title="seeds", y_axis="Sequenced Genome",
                 x_axis="Reference Genome", squeeze=False):
    plot = figure(title=title, plot_width=1000, plot_height=1000)

    if squeeze:
        x_squeezer = AxisSqueezer(squeeze)
        y_squeezer = AxisSqueezer(squeeze)
        pack_1 = Pack()
        pack_1.load(query_genome + "/ma/genome")
        pack_2 = Pack()
        pack_2.load(reference_genome + "/ma/genome")

        for idx in [*pack_1.contigStarts(), pack_1.unpacked_size_single_strand]:
            y_squeezer.add(idx)
        for idx in [*pack_2.contigStarts(), pack_2.unpacked_size_single_strand]:
            x_squeezer.add(idx)

        for y_start, seeds, _, _ in seeds_1:
            for seed in seeds:
                x_squeezer.add(seed.start_ref)
                x_squeezer.add(seed.start_ref + seed.size * (1 if seed.on_forward_strand else -1))
                y_squeezer.add(seed.start + y_start)
                y_squeezer.add(seed.start + seed.size + y_start)
        x_squeezer.prep()
        y_squeezer.prep()
        squeezed_seeds = []
        max_dist = 3
        for y_start__, seeds, a, b in seeds_1:
            for seed in combine_consecutive_diag_seeds(seeds):
                if seed.size == 0:
                    continue
                (y_start, y_end, x_start, x_end, forw) = (
                        y_squeezer.get(seed.start + y_start__),
                        y_squeezer.get(seed.start + seed.size + y_start__),
                        x_squeezer.get(seed.start_ref),
                        x_squeezer.get(seed.start_ref + seed.size * (1 if seed.on_forward_strand else -1)),
                        seed.on_forward_strand
                    )
                if len(squeezed_seeds) > 0 and squeezed_seeds[-1][4] == forw and abs(squeezed_seeds[-1][3] - x_start) <= max_dist:
                    squeezed_seeds[-1][1] = y_end
                    squeezed_seeds[-1][3] = x_end
                else:
                    squeezed_seeds.append([y_start, y_end, x_start, x_end, forw])

        x_squeezer_2 = AxisSqueezer(squeeze)
        y_squeezer_2 = AxisSqueezer(squeeze)
        for y_start, y_end, x_start, x_end, forw in squeezed_seeds:
            x_squeezer_2.add(x_start)
            x_squeezer_2.add(x_end)
            y_squeezer_2.add(y_start)
            y_squeezer_2.add(y_end)
        for idx in [*pack_1.contigStarts(), pack_1.unpacked_size_single_strand]:
            y_squeezer_2.add(y_squeezer.get(idx))
        for idx in [*pack_2.contigStarts(), pack_2.unpacked_size_single_strand]:
            x_squeezer_2.add(x_squeezer.get(idx))
        x_squeezer_2.prep()
        y_squeezer_2.prep()


        decorate_plot_helper(plot,
                            [x_squeezer_2.get(x_squeezer.get(i)) for i in [*pack_2.contigStarts(),
                                                                            pack_2.unpacked_size_single_strand]],
                            [*pack_2.contigNames()],
                            [y_squeezer_2.get(y_squeezer.get(i)) for i in [*pack_1.contigStarts(),
                                                                            pack_1.unpacked_size_single_strand]],
                            [*pack_1.contigNames()])

        squeezed_seeds_2 = [
            (
            y_squeezer_2.get(y_start),
            y_squeezer_2.get(y_end),
            x_squeezer_2.get(x_start),
            x_squeezer_2.get(x_end),
            forw)
            for y_start, y_end, x_start, x_end, forw in squeezed_seeds]
    else:
        decorate_plot(plot, query_genome, reference_genome)

    if True and not squeeze:
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
    if squeeze:
        for y_start, y_end, x_start, x_end, forw in squeezed_seeds_2:
            xs[False][forw].append([x_start, x_end])
            ys[False][forw].append([y_start, y_end])
            rs[False][forw].append(None)
    else:
        for y_start, seeds, _, filtered_seeds in seeds_1:
            for seedlist, filtered in [(seeds, False), (filtered_seeds, True)]:
                for seed_ in seedlist:
                    if filtered:
                        # uaaaagh @todo better construct needed
                        seed, reason = seed_
                    else:
                        seed = seed_
                        reason = None
                    if seed.size == 0:
                        continue
                    xs[filtered][seed.on_forward_strand].append([
                            seed.start_ref, seed.start_ref + seed.size * (1 if seed.on_forward_strand else -1)])
                    ys[filtered][seed.on_forward_strand].append([
                            seed.start + y_start, seed.start + seed.size + y_start])
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

genome_dir = main_data_folder + "/genome/yeasts/"

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

if True:
    out = []

    seeds_n_rects = compute_seeds(query_genome, reference_genome, db_name, seq_id)
    jumps = seeds_to_jumps(seeds_n_rects, query_genome, reference_genome)

    jumps_to_calls_to_db(jumps, db_name, query_genome, reference_genome)
    #exit()

    out.append(render_seeds(seeds_n_rects, query_genome, reference_genome))
    out.append(render_jumps(jumps, query_genome, reference_genome))
    show(row(out))

