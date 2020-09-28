from bokeh.plotting import figure, show, reset_output, ColumnDataSource
from bokeh.layouts import column, row
from bokeh.models import FuncTickFormatter
from bokeh.models.tools import HoverTool
from MA import *
from MSV import *
from sv_util.os_aligners import *
from bisect import bisect_right

global_prefix = "/MAdata/"

def load_genomes(query_genome, reference_genome):
    pack = Pack()
    pack.load(reference_genome + "/ma/genome")
    fm_index = FMIndex()
    #fm_index.load(reference_genome + "/ma/genome") # not used at the moment
    mm_index = MinimizerIndex(ParameterSetManager(), libMA.util.StringVector([str(pack.extract_forward_strand())]), libMA.util.StringVector(["chrVII"]))
    file_reader = FileReader(ParameterSetManager())
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
    param.by_name("Max Size Reseed").set(2000)
    param.by_name("Fixed SoC Width").set(100)
    param.by_name("Min NT in SoC").set(50)
    param.by_name("Rectangular SoC").set(False)

    pack, fm_index, mm_index, query_genomes = load_genomes(query_genome, reference_genome)
    mm_index.set_max_occ(ambiguity)

    seeding_module = MMFilteredSeeding(param)
    seed_lumper = SeedLumping(param)
    contig_filter = FilterContigBorder(param)
    soc_module = StripOfConsiderationSeeds(param)
    soc_filter = GetAllFeasibleSoCsAsSet(param)
    reseeding_1 = RecursiveReseedingSoCs(param, pack, 100)
    db_conn = DbConn({"SCHEMA": {"NAME": db_name}})
    mm_counter = HashFilterTable(db_conn).get_counter(seq_id)
    #mm_counter = HashCounter()
    jumps_from_seeds = SvJumpsFromSeeds(param, pack)
    reseeding = RecursiveReseeding(param, pack)

    ret = []
    for y_start, query_genome in query_genomes:
        seeds = seeding_module.execute(mm_index, query_genome, pack, mm_counter)
        mems = seed_lumper.execute(seeds, query_genome, pack)
        #ret.append((y_start, mems, []))
        #continue
        filtered_mems = contig_filter.execute(mems, pack)
        socs = soc_module.execute(filtered_mems, query_genome, pack)
        filtered_seeds = soc_filter.execute(socs)
        filtered_seeds_2 = reseeding_1.execute(filtered_seeds, pack, query_genome)
        helper_ret = jumps_from_seeds.cpp_module.execute_helper(filtered_seeds_2, pack, query_genome)
        #reseeded_mems = helper_ret.seeds
        #layer_of_seeds = helper_ret.layer_of_seeds
        rectangles = helper_ret.rectangles
        parlindromes = helper_ret.parlindrome
        #fill_of_rectangles = helper_ret.rectangles_fill
        #seed_sample_sizes = helper_ret.rectangle_ambiguity
        #rectangle_used_dp = helper_ret.rectangle_used_dp

        #ret.append((y_start, filtered_seeds_2, rectangles))
        ret.append((y_start, reseeding.execute(filtered_seeds_2, pack, query_genome), rectangles))
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
        seeds_list = []
        while not f_stream.eof():
            alignment = file_reader.execute(f_stream, pack, read)
            seeds = alignment.to_seeds(pack)
            for seed in seeds:
                seeds_list.append(seed)

        ret.append((y_start, seeds_list))

    return ret


class SeedsFilter:
    def __init__(self, intervals=None, intervals_2=None):
        self.intervals = intervals
        if not self.intervals is None:
            self.keys = [start for start, end in intervals]
        self.intervals_2 = intervals_2
        if not self.intervals_2 is None:
            self.keys_2 = [start for start, end in intervals_2]

    def __filter(self, s_start, s_end, keys, intervals):
        idx = bisect_right(keys, s_start)
        if idx != 0:
            start, end = intervals[idx-1]
            s = max(start, s_start)
            e = min(end, s_end)
            if e - s > 0:
                return True
        if idx < len(intervals):
            start, end = intervals[idx]
            s = max(start, s_start)
            e = min(end, s_end)
            if e - s > 0:
                return True
        return False

    def seed_filtered(self, seed):
        if seed.on_forward_strand:
            s_start = seed.start_ref
            s_end = seed.start_ref + seed.size
        else:
            s_start = seed.start_ref - seed.size
            s_end = seed.start_ref
        ret = False
        if not self.intervals is None:
            ret = ret or self.__filter(s_start, s_end, self.keys, self.intervals)
        if not self.intervals_2 is None:
            ret = ret or self.__filter(seed.start, seed.start + seed.size, self.keys_2, self.intervals_2)
        return ret

def render_seeds(seeds, merged_intervals=None, merged_intervals_2=None):
    plot = figure(title="seeds", plot_width=1000, plot_height=1000)
    xs = {True:{True:[], False:[]}, False:{True:[], False:[]}}
    ys = {True:{True:[], False:[]}, False:{True:[], False:[]}}
    cs = {True:{True:"lightblue", False:"yellow"}, False:{True:"blue", False:"orange"}}
    if not merged_intervals is None:
        seed_filter = SeedsFilter(merged_intervals, merged_intervals_2)
        ref_size = max([seed.start + seed.size for seed in seeds])
        s = [start for start, end in merged_intervals]
        e = [end for start, end in merged_intervals]
        plot.quad(left="start", right="end", top=ref_size, bottom=0, fill_color="pink", line_width=0, 
                  source=ColumnDataSource(data={ "start":s, "end":e }))
    if not merged_intervals_2 is None:
        ref_size = max([seed.start_ref + seed.size for seed in seeds])
        s = [start for start, end in merged_intervals_2]
        e = [end for start, end in merged_intervals_2]
        plot.quad(bottom="start", top="end", right=ref_size, left=0, fill_color="pink", line_width=0, 
                  source=ColumnDataSource(data={ "start":s, "end":e }))
    for seed in seeds:
        filtered = False
        if not merged_intervals is None:
            filtered = seed_filter.seed_filtered(seed)
        xs[filtered][seed.on_forward_strand].append(seed.start_ref)
        xs[filtered][seed.on_forward_strand].append(seed.start_ref + seed.size * (1 if seed.on_forward_strand else -1))
        ys[filtered][seed.on_forward_strand].append(seed.start)
        ys[filtered][seed.on_forward_strand].append(seed.start + seed.size)
        xs[filtered][seed.on_forward_strand].append(float("NaN"))
        ys[filtered][seed.on_forward_strand].append(float("NaN"))
    for filtered in [False]:
        for forw in [True, False]:
            plot.line(x="xs", y="ys", color=cs[filtered][forw], line_width=4, line_cap="round", source=ColumnDataSource(data={
                    "xs":xs[filtered][forw], "ys":ys[filtered][forw]
                }))
    return plot

def seeds_to_jumps(seeds_n_rects, query_genome, reference_genome):
    param = ParameterSetManager()
    pack, fm_index, mm_index, query_genomes = load_genomes(query_genome, reference_genome)
    jumps_from_seeds = SvJumpsFromExtractedSeeds(param, pack)
    filter_by_ambiguity = FilterJumpsByRefAmbiguity(param)

    ret = []
    for (_, seeds, _), (_, query_genome) in zip(seeds_n_rects, query_genomes):
        jumps = jumps_from_seeds.execute(libMA.containers.Seeds(seeds), pack, query_genome)
        filtered_jumps = filter_by_ambiguity.execute(jumps, pack)
        ret.append(jumps)
    return ret

def decorate_plot(plot, query_genome, reference_genome, diagonal=False):
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
        ys.append(pack_1.unpacked_size_single_strand)
        xs.append(float("NaN"))
        ys.append(float("NaN"))
    plot.line(x=xs, y=ys, color="black", line_width=1)
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
    idx = 0
    colors = {True:{True:"blue", False:"green"}, False:{True:"purple", False:"orange"}}
    for jumps_ in jumps:
        for jump in jumps_:
            #if idx != 8281:
            #    idx+=1
            #    continue
            f = jump.from_pos
            t = jump.to_pos
            if not jump.from_known():
                f = t
                t += 1
            if not jump.to_known():
                t = f
                f -= 1
            if abs(f - t) < 200:
                continue
            xs.append(f)
            ys.append(t)
            ms.append(jump.was_mirrored)
            if jump.query_to - jump.query_from >= 5000:
                cs.append("red")
            else:
                cs.append(colors[jump.from_forward][jump.to_forward])
            ids.append(idx)
            qlens.append(jump.query_to - jump.query_from )
            fs.append("fow" if jump.from_forward else "rev")
            ts.append("fow" if jump.to_forward else "rev")
            idx+=1
    plot = figure(title="entries", plot_width=1000, plot_height=1000)
    decorate_plot(plot, reference_genome, reference_genome, True)
    plot.x(x="xs", y="ys", color="cs", line_width=4, source=ColumnDataSource(data={
                    "xs":xs, "ys":ys, "cs":cs, "ms":ms, "ids": ids, "qlens": qlens, "fs":fs, "ts":ts
                }))
    plot.add_tools(HoverTool(tooltips=[("x, y", "@xs, @ys"),("id", "@ids"),("color, mirrored", "@cs, @ms"),
                                       ("query dist", "@qlens"), ("strandinfo", "from @fs to @ts")]))
    return plot

def jumps_to_calls_to_db(jumps, db_name, query_genome_str, reference_genome, max_q_dist=None, min_dist=None):
    db_conn = DbConn({"SCHEMA": {"NAME": db_name}})

    JumpRunTable(db_conn)
    SvCallerRunTable(db_conn)
    name = "Ground Truth"
    if not max_q_dist is None:
        name = name + " Q-Dist-Filter:" + str(max_q_dist)
    get_inserter = GetCallInserter(ParameterSetManager(), db_conn, name,
                                   "SV's form genome assembly", -1)
    pool = PoolContainer(1, db_name)
    sv_inserter = get_inserter.execute(pool)
    
    pack, fm_index, mm_index, query_genomes = load_genomes(query_genome_str, reference_genome)

    def at_contig_border(x, dist=3):
        idx = pack.seq_id_for_pos(x)
        if pack.start_of_sequence_id(idx) + dist >= x:
            return True
        if pack.end_of_sequence_id(idx) <= x + dist:
            return True
        return False

    idx = 0
    for jumps_, (_, query_genome) in zip(jumps, query_genomes):
        query_genome.check()
        for jump in jumps_:
            f = jump.from_pos
            t = jump.to_pos
            if not jump.from_known():
                continue
            if not jump.to_known():
                continue
            if at_contig_border(f) or at_contig_border(t):
                continue
            if not max_q_dist is None:
                if jump.query_to - jump.query_from > max_q_dist:
                    continue 
            if not min_dist is None:
                if jump.query_to - jump.query_from < min_dist and abs(t - f) < min_dist:
                    continue
            call = SvCall(f, t, 0, 0, jump.from_forward, jump.to_forward, 1000)
            if jump.query_from < jump.query_to:
                call.inserted_sequence = NucSeq(query_genome, jump.query_from, jump.query_to)
                call.inserted_sequence.check()
            call.order_id = idx
            idx += 1
            call.mirrored = jump.was_mirrored
            sv_inserter.insert(call)

    sv_inserter.close(pool)

    return get_inserter.cpp_module.id, None

def render_seeds_2(seeds_1, seeds_2, query_genome, reference_genome, title="seeds"):
    plot = figure(title=title, plot_width=1000, plot_height=1000)

    decorate_plot(plot, query_genome, reference_genome)

    xs = []
    xe = []
    ys = []
    ye = []
    for y_start, _, rects in seeds_1:
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
    cs = {True:{True:"purple", False:"green"}, False:{True:"blue", False:"orange"}}
    lw = {True:8, False: 4}
    for y_start, seeds, _ in seeds_1:
        for seed in seeds:
            if seed.size < 20:
                continue
            filtered = False
            xs[filtered][seed.on_forward_strand].append(seed.start_ref)
            xs[filtered][seed.on_forward_strand].append(seed.start_ref + seed.size * (1 if seed.on_forward_strand else -1))
            ys[filtered][seed.on_forward_strand].append(seed.start + y_start)
            ys[filtered][seed.on_forward_strand].append(seed.start + seed.size + y_start)
            xs[filtered][seed.on_forward_strand].append(float("NaN"))
            ys[filtered][seed.on_forward_strand].append(float("NaN"))
    if not seeds_2 is None:
        for y_start, seeds, _ in seeds_2:
            for seed in seeds:
                filtered = True
                xs[filtered][seed.on_forward_strand].append(seed.start_ref)
                xs[filtered][seed.on_forward_strand].append(seed.start_ref + seed.size * (1 if seed.on_forward_strand else -1))
                ys[filtered][seed.on_forward_strand].append(seed.start + y_start)
                ys[filtered][seed.on_forward_strand].append(seed.start + seed.size + y_start)
                xs[filtered][seed.on_forward_strand].append(float("NaN"))
                ys[filtered][seed.on_forward_strand].append(float("NaN"))
    for filtered in [True, False]:
        for forw in [True, False]:
            plot.line(x="xs", y="ys", color=cs[filtered][forw], line_width=lw[filtered], line_cap="round", source=ColumnDataSource(data={
                    "xs":xs[filtered][forw], "ys":ys[filtered][forw]
                }))

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

if __name__ == "__main__":

    seeds_n_rects = compute_seeds(query_genome, reference_genome, "UFRJ50816_test_reconstruct", 1)
    jumps = seeds_to_jumps(seeds_n_rects, query_genome, reference_genome)
    #jumps_to_calls_to_db(jumps, "UFRJ50816_test_reconstruct", query_genome, reference_genome)
    #jumps_to_calls_to_db(jumps, "UFRJ50816_test_reconstruct", query_genome, reference_genome, 5000)
    exit()

    aligner_seeds = None
    #aligner_seeds = run_aligner(query_genome, reference_genome)
    #aligner_seeds = [(x,y,[]) for x,y in aligner_seeds]

    out = []
    out.append(render_seeds_2(seeds_n_rects, aligner_seeds, query_genome, reference_genome))
    out.append(render_jumps(jumps, query_genome, reference_genome))
    show(row(out))


# ~/workspace/samtools/samtools view -h minimap2.sorted.bam "chrIX:251682-259682" > minimap2.filtered.bam


# ~/workspace/samtools/samtools view minimap2.filtered.bam | awk '{print length($10), $1}' | sort -n


# "chrVII:778311-788311"

# m150417_101358_00127_c100782972550000001823173608251533_s1_p0/39217/0_22799

#1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/44963/0_5594;1:m150423_092113_00127_c100802542550000001823174910081570_s1_p0/109922/8352_13954;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/18037/0_5603;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/18037/5658_11267;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/20956/20700_26370;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/144875/5623_11295;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/71304/464_11950;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/13415/20843_26706;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/36188/5946_11812;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/36188/0_5903;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/36188/11858_17775;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/36188/17820_23740;1:m150423_092113_00127_c100802542550000001823174910081570_s1_p0/141524/5704_11651;1:m150423_092113_00127_c100802542550000001823174910081570_s1_p0/69336/6856_12806;1:m150423_092113_00127_c100802542550000001823174910081570_s1_p0/53772/16676_22651;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/48812/0_6104;

# 1:m150423_092113_00127_c100802542550000001823174910081570_s1_p0/69336/579_6809;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/142877/1520_8074;1:m150424_132538_00127_c100802542550000001823174910081577_s1_p0/56174/0_6597;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/144999/4834_11460;1:m150424_132538_00127_c100802542550000001823174910081577_s1_p0/96953/0_6679;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/119223/0_6698;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/89902/0_6727;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/84553/0_6764;1:m150424_132538_00127_c100802542550000001823174910081577_s1_p0/96953/6725_13515;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/13415/6922_13772;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/13415/0_6878;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/13415/13818_20797;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/67617/0_14819;1:m150424_132538_00127_c100802542550000001823174910081577_s1_p0/96953/20946_28172;1:m150424_132538_00127_c100802542550000001823174910081577_s1_p0/131853/3711_10942;1:m150424_132538_00127_c100802542550000001823174910081577_s1_p0/96953/13559_20900;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/52958/0_7342;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/114215/19087_26459;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/153159/23675_31070;

# 1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/50652/0_7727;1:m150424_132538_00127_c100802542550000001823174910081577_s1_p0/115600/2203_10049;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/144875/11342_19397;1:m150423_092113_00127_c100802542550000001823174910081570_s1_p0/23917/0_8142;1:m150424_132538_00127_c100802542550000001823174910081577_s1_p0/107598/0_8176;1:m150423_092113_00127_c100802542550000001823174910081570_s1_p0/53772/0_8260;1:m150423_092113_00127_c100802542550000001823174910081570_s1_p0/109922/0_8306;1:m150423_092113_00127_c100802542550000001823174910081570_s1_p0/53772/8308_16631;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/97715/0_8460;1:m150423_092113_00127_c100802542550000001823174910081570_s1_p0/142063/326_8861;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/2179/0_8554;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/2179/8599_17162;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/140135/911_9670;1:m150423_092113_00127_c100802542550000001823174910081570_s1_p0/134366/0_8891;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/98776/9081_18050;1:m150423_092113_00127_c100802542550000001823174910081570_s1_p0/18417/786_9783;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/98776/0_9035;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/126913/19910_29096;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/108148/0_9205;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/146893/11756_21019;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/52204/0_9265;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/114215/0_9279;1:m150424_132538_00127_c100802542550000001823174910081577_s1_p0/109222/0_9291;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/145770/5355_14713;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/122350/0_9424;1:m150423_092113_00127_c100802542550000001823174910081570_s1_p0/145964/0_9575;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/126913/10149_19866;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/114215/9322_19044;1:m150423_092113_00127_c100802542550000001823174910081570_s1_p0/47660/3687_13459;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/61988/0_9778;1:m150424_132538_00127_c100802542550000001823174910081577_s1_p0/161131/21667_31448;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/157114/5062_14887;1:m150424_132538_00127_c100802542550000001823174910081577_s1_p0/119359/3182_13040;1:m150424_132538_00127_c100802542550000001823174910081577_s1_p0/130318/0_9884;1:m150424_132538_00127_c100802542550000001823174910081577_s1_p0/49524/0_10347;1:m150423_092113_00127_c100802542550000001823174910081570_s1_p0/386/0_10669;1:m150424_132538_00127_c100802542550000001823174910081577_s1_p0/59588/0_11162;1:m150423_092113_00127_c100802542550000001823174910081570_s1_p0/55563/11668_23254;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/66728/13702_25321;1:m150423_092113_00127_c100802542550000001823174910081570_s1_p0/55563/0_11621;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/95246/0_11806;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/129268/0_11829;1:m150423_092113_00127_c100802542550000001823174910081570_s1_p0/139555/0_12158;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/150846/0_12405;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/9765/0_12518;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/105516/0_13072;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/129509/0_13118;1:m150424_132538_00127_c100802542550000001823174910081577_s1_p0/105081/13820_27418;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/66728/0_13658;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/117015/0_13993;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/131549/0_14564;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/67617/0_14819;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/134928/0_15752;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/35421/2824_18781;1:m150423_092113_00127_c100802542550000001823174910081570_s1_p0/113604/0_16179;1:m150424_132538_00127_c100802542550000001823174910081577_s1_p0/155075/0_17159;1:m150424_132538_00127_c100802542550000001823174910081577_s1_p0/127492/0_19215;1:m150423_092113_00127_c100802542550000001823174910081570_s1_p0/146511/1919_21890;1:m150424_132538_00127_c100802542550000001823174910081577_s1_p0/161991/0_20332;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/20956/0_20653;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/39217/0_22799;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/50073/0_23026


# "chrXV:1065771-1077771"

#1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/55885/0_8619;1:m150424_132538_00127_c100802542550000001823174910081577_s1_p0/35295/0_8639;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/129742/17067_25729;1:m150423_092113_00127_c100802542550000001823174910081570_s1_p0/152119/0_8842;1:m150424_132538_00127_c100802542550000001823174910081577_s1_p0/114840/7346_16216;1:m150424_132538_00127_c100802542550000001823174910081577_s1_p0/109522/13296_22292;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/105542/33101_42099;1:m150424_132538_00127_c100802542550000001823174910081577_s1_p0/70172/11458_20725;1:m150423_092113_00127_c100802542550000001823174910081570_s1_p0/102214/9526_18885;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/41024/0_9412;1:m150424_132538_00127_c100802542550000001823174910081577_s1_p0/70172/1965_11411;1:m150423_092113_00127_c100802542550000001823174910081570_s1_p0/102214/0_9481;1:m150423_092113_00127_c100802542550000001823174910081570_s1_p0/135212/0_9492;1:m150423_092113_00127_c100802542550000001823174910081570_s1_p0/30757/0_9629;1:m150424_132538_00127_c100802542550000001823174910081577_s1_p0/154519/0_9866;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/10228/6379_16379;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/33531/10262_20295;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/17978/10250_20357;1:m150423_092113_00127_c100802542550000001823174910081570_s1_p0/12964/2375_12542;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/17978/0_10208;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/33531/0_10215;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/74182/0_10309;1:m150424_132538_00127_c100802542550000001823174910081577_s1_p0/100090/0_10419;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/135130/0_10422;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/56647/2530_12954;1:m150424_132538_00127_c100802542550000001823174910081577_s1_p0/126768/0_10441;1:m150424_132538_00127_c100802542550000001823174910081577_s1_p0/100421/0_10472;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/141610/15767_26292;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/126483/15105_25714;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/56864/0_10624;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/56864/10667_21319;1:m150424_132538_00127_c100802542550000001823174910081577_s1_p0/100421/10515_21179;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/52520/14743_25415;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/89986/0_10952;1:m150423_092113_00127_c100802542550000001823174910081570_s1_p0/105751/0_11049;1:m150423_092113_00127_c100802542550000001823174910081570_s1_p0/39529/2192_13268;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/138754/0_11288;1:m150423_092113_00127_c100802542550000001823174910081570_s1_p0/26078/0_11407;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/123022/0_11429;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/37548/17257_28842;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/83276/0_11654;1:m150423_092113_00127_c100802542550000001823174910081570_s1_p0/41545/12928_24890;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/126440/0_12051;1:m150423_092113_00127_c100802542550000001823174910081570_s1_p0/12986/0_12256;1:m150424_132538_00127_c100802542550000001823174910081577_s1_p0/137188/0_12376;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/130987/8198_20577;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/121636/0_12585;1:m150423_092113_00127_c100802542550000001823174910081570_s1_p0/41545/0_12885;1:m150423_092113_00127_c100802542550000001823174910081570_s1_p0/3142/0_13165;1:m150423_092113_00127_c100802542550000001823174910081570_s1_p0/26833/13390_26583;1:m150424_132538_00127_c100802542550000001823174910081577_s1_p0/109522/0_13253;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/84441/1782_15048;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/141043/13399_26728;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/141043/0_13354;1:m150424_132538_00127_c100802542550000001823174910081577_s1_p0/16006/7084_20531;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/56436/8436_22097;1:m150423_092113_00127_c100802542550000001823174910081570_s1_p0/87463/0_14608;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/52520/0_14698;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/144355/1327_16770;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/156966/0_15531;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/141610/0_15722;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/15467/0_16528;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/116614/0_16703;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/49306/0_16755;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/129763/11015_27846;1:m150424_132538_00127_c100802542550000001823174910081577_s1_p0/95551/0_17024;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/37548/0_17210;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/95743/0_17357;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/20949/0_18104;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/135400/0_18733;1:m150423_092113_00127_c100802542550000001823174910081570_s1_p0/52061/11034_30036;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/32549/0_20061;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/50125/0_25490;1:m150416_014538_00127_c100782992550000001823173608251512_s1_p0/94848/0_26413;1:m150417_101358_00127_c100782972550000001823173608251533_s1_p0/105488/0_27221;


# m150417_101358_00127_c100782972550000001823173608251533_s1_p0/36119/0_25862