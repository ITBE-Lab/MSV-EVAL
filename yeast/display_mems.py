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
    param.set_selected("SV-PacBio")

    pack, fm_index, mm_index, query_genomes = load_genomes(query_genome, reference_genome)
    mm_index.set_max_occ(ambiguity)

    seeding_module = MMFilteredSeeding(param)
    seed_lumper = SeedLumping(param)
    contig_filter = FilterContigBorder(param)
    soc_module = StripOfConsiderationSeeds(param)
    soc_filter = GetAllFeasibleSoCsAsSet(param)
    reseeding_1 = RecursiveReseedingSoCs(param, pack)
    db_conn = DbConn({"SCHEMA": {"NAME": db_name}})
    mm_counter = HashFilterTable(db_conn).get_counter(seq_id)
    #mm_counter = HashCounter()
    jumps_from_seeds = SvJumpsFromSeeds(param, pack)
    #reseeding = RecursiveReseeding(param, pack)

    ret = []
    for y_start, query_genome in query_genomes:
        seeds = seeding_module.execute(mm_index, query_genome, pack, mm_counter)
        mems = seed_lumper.execute(seeds, query_genome, pack)
        #ret.append((y_start, mems, []))
        #continue
        filtered_mems = contig_filter.execute(mems, pack)
        socs = soc_module.execute(filtered_mems, query_genome, pack)
        soc_filtered_seeds = soc_filter.execute(socs)
        filtered_seeds_2, helper_ret_1 = reseeding_1.cpp_module.execute_helper(soc_filtered_seeds, pack, query_genome)
        #helper_ret, seeds_2 = jumps_from_seeds.cpp_module.execute_helper(filtered_seeds_2, pack, query_genome)
        #reseeded_mems = helper_ret.seeds
        #layer_of_seeds = helper_ret.layer_of_seeds
        #rectangles = helper_ret_1.rectangles
        #parlindromes = helper_ret.parlindrome
        #overlapping = helper_ret.overlapping
        #fill_of_rectangles = helper_ret.rectangles_fill
        #seed_sample_sizes = helper_ret.rectangle_ambiguity
        #rectangle_used_dp = helper_ret.rectangle_used_dp

        filtered_seeds = []
        while not socs.empty():
            for seed in socs.pop():
                filtered_seeds.append( (seed, "initial SoC") )
        for seed in helper_ret_1.seed_removed:
            filtered_seeds.append( (seed, "reseeding Soc / overlapping filter / enclosed SoC") )
        for seed, parlindrome, overlapping in zip(helper_ret_1.seeds, helper_ret_1.parlindrome,
                                                  helper_ret_1.overlapping):
            if parlindrome:
                filtered_seeds.append((seed, "palrindrome"))
            if overlapping:
                filtered_seeds.append((seed, "overlapping"))


        #ret.append((y_start, filtered_seeds_2, rectangles))
        ret.append((y_start, filtered_seeds_2, helper_ret_1.rectangles, filtered_seeds))
        #ret.append((y_start, reseeding.execute(filtered_seeds_2, pack, query_genome), rectangles))
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

def make_read_range_table(db_name, seq_id, assembled_genome, use_mm2=True):
    read_set = {"technology":"pb", "name":"test", "fasta_file":"reads.fasta"}
    json_dict = {"reference_path":reference_genome}
    sam_file_path = "read_ext_out.sam"
    ret = []
    with open("reads.fasta", "w") as out_file:
        for read in iterate_reads(ParameterSetManager(), db_name, seq_id):
            out_file.write(">")
            out_file.write(str(read.id))
            out_file.write("\n")
            out_file.write(str(read))
            out_file.write("\n")

    if use_mm2:
        mm2(read_set, sam_file_path, json_dict)
    else:
        bwa_single(read_set, sam_file_path, json_dict)
    f_stream = FileStreamFromPath("read_ext_out.sam")
    read_by_name = ReadByName()
    # just always return nullptr since we do not need the cigars
    read_by_name.return_null_for_unknown = True
    file_reader = SamFileReader(ParameterSetManager())
    db_conn = DbConn({"SCHEMA": {"NAME": db_name}})
    ext_table = ReadRangeTable(db_conn)
    ext_table.drop_indices()
    pack = Pack()
    pack.load(assembled_genome + "/ma/genome")
    while not f_stream.eof():
        alignment = file_reader.execute(f_stream, pack, read_by_name)
        read_id = int(alignment.stats.name)
        ext_table.insert(read_id, alignment)
    ext_table.gen_indices()

def make_read_range_table_from_simulated_n_seeds(db_name, seq_id, assembled_genome, reference_genome):
    db_conn = DbConn({"SCHEMA": {"NAME": db_name}})
    ext_table = ReadRangeTable(db_conn)
    ext_table.drop_indices()
    pack = Pack()
    pack.load(assembled_genome + "/ma/genome")
    pack_2 = Pack()
    pack_2.load(reference_genome + "/ma/genome")

    parameter_set_manager = ParameterSetManager()
    parameter_set_manager.set_selected("SV-PacBio")
    mm_index = MinimizerIndex(parameter_set_manager, pack.contigSeqs(), pack.contigNames())
    mm_index.set_max_occ(2)
    mm_counter = HashFilterTable(db_conn).get_counter(seq_id)
    seeding_module = MMFilteredSeeding(parameter_set_manager)
    seed_lumper = SeedLumping(parameter_set_manager)
    contig_filter = FilterContigBorder(parameter_set_manager)
    soc_module = StripOfConsiderationSeeds(parameter_set_manager)
    soc_filter = GetAllFeasibleSoCsAsSet(parameter_set_manager)
    read_table = ReadTable(db_conn)

    with open("reads.fasta", "w") as out_file:
        for read in iterate_reads(ParameterSetManager(), db_name, seq_id):
            split = read_table.read_name(read.id).split("_")
            #print(split)
            chrom = split[0]
            start_pos = int(split[1]) + pack.start_of_sequence("chr" + chrom[3:])
            strand = split[2]
            end_pos = start_pos + len(read)
            ext_table.insert_range(read, start_pos, end_pos, 1, False)

            minimizers = seeding_module.execute(mm_index, read, pack, mm_counter)
            lumped_seeds = seed_lumper.execute(minimizers, read, pack)
            c_filter_seeds = contig_filter.execute(lumped_seeds, pack)
            socs = soc_module.execute(c_filter_seeds, read, pack)
            filtered_seeds_pledge_2 = soc_filter.execute(socs)
            start_pos_mem = end_pos
            end_pos_mem = start_pos

            for mems in filtered_seeds_pledge_2.content:
                for mem in mems:
                    mem_start = mem.start_ref
                    mem_end = mem.start_ref + mem.size
                    if not mem.on_forward_strand:
                        mem_start = mem.start_ref - mem.size
                        mem_end = mem.start_ref

                    if mem_start < end_pos and mem_end > start_pos:
                        start_pos_mem = min(start_pos_mem, mem_start)
                        end_pos_mem = max(end_pos_mem, mem_end)

            if start_pos_mem < end_pos_mem:
                ext_table.insert_range(read, start_pos_mem, end_pos_mem, 1, True)

    ext_table.gen_indices()

def view_coverage(db_name, seq_id, assembled_genome, steps=10000, with_selection=False):
    db_conn = DbConn({"SCHEMA": {"NAME": db_name}})
    ext_table = ReadRangeTable(db_conn, with_selection)
    ext_table.gen_indices() # make sure indices exist (does nothing if they are present)
    pack = Pack()
    pack.load(assembled_genome + "/ma/genome")
    plot = figure(title="coverage of " + assembled_genome, plot_width=1000, plot_height=1000)
    l = [
        (False, None, "simulated reads", "red"),
        #(True, 0, "primary & map_q>=0", "blue"),
        (True, None, "unambiguous simulated reads", "green")
    ]
    max_y = 0
    for primary, map_q_min, name, color in l:
        xs = [0]
        ys = [0]
        for x in range(0, pack.unpacked_size_single_strand, max(1, pack.unpacked_size_single_strand // steps)):
            xs.append(x)
            if map_q_min is None:
                y = ext_table.coverage(x, x+1, seq_id, primary)
            else:
                y = ext_table.coverage(x, x+1, seq_id, primary, map_q_min)
            max_y = max(y, max_y)
            ys.append(y)
        xs.append(pack.unpacked_size_single_strand)
        ys.append(0)
        plot.patch(x=xs, y=ys, legend_label=name, color=color)

    #for x in [*pack.contigStarts(), pack.unpacked_size_single_strand]:
    #    plot.line(x=[x, x], y=[0, max_y], color="black")

    decorate_plot(plot, None, assembled_genome, x_only=True, max_y=max_y*1.25)
    plot.xaxis.axis_label = "Sequenced Genome"
    plot.yaxis.axis_label = "Coverage"

    return plot

def reads_per_jump(db_name, jumps, seq_id, min_stick_out=50, primary=True, min_map_q=0, min_cov=10):
    db_conn = DbConn({"SCHEMA": {"NAME": db_name}})
    ext_table = ReadRangeTable(db_conn)
    ext_table.gen_indices() # make sure indices exist (does nothing if they are present)

    ret = []
    for idx, jumps_ in enumerate(jumps):
        print("reads_per_jump:", idx, "/", len(jumps), "...")
        ret.append([])
        for jump in jumps_:
            if min_map_q is None:
                cov = ext_table.coverage(jump.query_from - min_stick_out, jump.query_to + min_stick_out, seq_id,
                                            primary)
            else:
                cov = ext_table.coverage(jump.query_from - min_stick_out, jump.query_to + min_stick_out, seq_id,
                                            primary, min_map_q)
            if cov < min_cov or True:
                if min_map_q is None:
                    cov_start = ext_table.coverage(jump.query_from - min_stick_out, jump.query_from + min_stick_out,
                                                    seq_id, primary)
                    cov_end = ext_table.coverage(jump.query_to - min_stick_out, jump.query_to + min_stick_out, seq_id,
                                                primary)
                else:
                    cov_start = ext_table.coverage(jump.query_from - min_stick_out, jump.query_from + min_stick_out,
                                                    seq_id, primary, min_map_q)
                    cov_end = ext_table.coverage(jump.query_to - min_stick_out, jump.query_to + min_stick_out, seq_id,
                                                primary, min_map_q)
            else:
                cov_start = None
                cov_end = None
            ret[-1].append( (cov, cov_start, cov_end) )
    return ret

def render_reads_per_jump(jumps, cov_lists, assembled_genome):
    plot = figure(title="#reads that cover a call " + assembled_genome, plot_width=1000, plot_height=1000)
    plot4 = figure(title="#reads that cover a breakpoint " + assembled_genome, plot_width=1000, plot_height=1000)
    plot2 = figure(title="query distance of calls (positive: cov < 10; negative: cov >= 10) of " + assembled_genome, plot_width=1000, plot_height=1000)

    def dist_list_to_bars(dist_list, plot, name, color, positive=True, num_bars=100):
        dist_list.sort()
        buckets = [0]*(num_bars+1)
        bucket_len = (dist_list[-1] - dist_list[0]) / num_bars
        for dist in dist_list:
            buckets[ int((dist - dist_list[0]) / bucket_len) ] += 1 if positive else -1
        xs = [x * bucket_len + bucket_len/2 for x in range(num_bars+1)]
        plot.vbar(x=xs, bottom=0 if positive else buckets, top=buckets if positive else 0,
                  width=bucket_len*0.75, legend_label=name, color=color, fill_alpha=.5)

    def cov_patch(cov_dict, plot, name, color):
        xs = [0]
        ys = [0]
        #print("coverage_dict ( coverage -> (num,id) ):")
        for (cov, l) in sorted(list(cov_dict.items())):
            #print(cov, len(l), l)
            xs.append(cov)
            ys.append(len(l))

        xs.append(cov)
        ys.append(0)

        plot.patch(x=xs, y=ys, legend_label=name, color=color, fill_alpha=.5)

    for (cov_list, name, color) in cov_lists:
        cov_dict = {}
        cov_dict_bp = {}
        idx = 0
        for cov_list_, jumps_ in zip(cov_list, jumps):
            for (cov, cov_start, cov_end), jump in zip(cov_list_, jumps_):
                if not cov in cov_dict:
                    cov_dict[cov] = []
                cov_dict[cov].append( (idx, jump.query_to - jump.query_from) )

                if not cov_start in cov_dict_bp:
                    cov_dict_bp[cov_start] = []
                cov_dict_bp[cov_start].append(idx)

                if not cov_end in cov_dict_bp:
                    cov_dict_bp[cov_end] = []
                cov_dict_bp[cov_end].append(idx)

                idx += 1

        cov_patch(cov_dict, plot, name, color)
        cov_patch(cov_dict_bp, plot4, name, color)

        q_dist_list = []
        q_dist_list2 = []
        for (cov, l) in sorted(list(cov_dict.items())):
            for idx, dist in l:
                if cov < 10:
                    q_dist_list.append(dist)
                else:
                    q_dist_list2.append(dist)
        dist_list_to_bars(q_dist_list, plot2, name, color)
        dist_list_to_bars(q_dist_list2, plot2, name, color, False)

    plot.xaxis.axis_label = "coverage"
    plot.yaxis.axis_label = "#calls"

    plot2.xaxis.axis_label = "q_dist"
    plot2.yaxis.axis_label = "#calls"


    return [plot, plot4, plot2]

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

def render_jumps(jumps, jump_cov, query_genome, reference_genome, min_cov=10):
    xs = []
    ys = []
    cs = []
    ms = []
    ids = []
    covs = []
    cov_starts = []
    cov_ends = []
    qlens = []
    fs = []
    ts = []
    colors = {True:{True:"blue", False:"green"}, False:{True:"purple", False:"orange"}}
    for jumps_, jump_cov_ in zip(jumps, jump_cov):
        for jump, (cov, cov_start, cov_end) in zip(jumps_, jump_cov_):
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
            #if abs(f - t) < 200:
            #    continue
            xs.append(f)
            ys.append(t)
            ms.append(jump.was_mirrored)
            covs.append(cov)
            cov_starts.append(cov_start)
            cov_ends.append(cov_end)
            if cov < min_cov:
                if cov_start >= min_cov and cov_end >= min_cov:
                    cs.append("black")
                else:
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
                    "xs":xs, "ys":ys, "cs":cs, "ms":ms, "ids": ids, "qlens": qlens, "fs":fs, "ts":ts, "cov": covs,
                    "cov_start": cov_starts, "cov_end":cov_ends
                }))
    plot.add_tools(HoverTool(tooltips=[("x, y", "@xs, @ys"),("id", "@ids"),("color, mirrored", "@cs, @ms"),
                                       ("query dist", "@qlens"), ("strandinfo", "from @fs to @ts"),
                                       ("read coverage", "enclose: @cov, start: @cov_start, end: @cov_end")]))
    plot.xaxis.axis_label = "Reference Genome"
    plot.yaxis.axis_label = "Reference Genome"
    return plot

def filter_jumps(jumps, reference_genome, max_q_dist=None, min_dist=None):
    pack = Pack()
    pack.load(reference_genome + "/ma/genome")
    def at_contig_border(x, dist=3000):
        idx = pack.seq_id_for_pos(x)
        if pack.start_of_sequence_id(idx) + dist >= x:
            return True
        if pack.end_of_sequence_id(idx) <= x + dist:
            return True
        return False

    ret = []
    for jumps_ in jumps:
        ret.append([])
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
            ret[-1].append(jump)
    return ret

def jumps_to_calls_to_db(jumps, cov_list, db_name, query_genome_str, reference_genome, min_cov=10):
    parameter_set_manager = ParameterSetManager()
    parameter_set_manager.set_selected("SV-PacBio")
    min_size = parameter_set_manager.by_name("Min Size Edge").get()
    db_conn = DbConn({"SCHEMA": {"NAME": db_name}})

    JumpRunTable(db_conn)
    SvCallerRunTable(db_conn)
    sv_call_table = SvCallTable(db_conn)
    one_sided_calls = OneSidedCallsTable(db_conn)
    caller_run_id = GetCallInserter(parameter_set_manager, db_conn, "Ground Truth",
                                   "SV's form genome assembly", -1).cpp_module.id
    low_cov_caller_run_id = GetCallInserter(parameter_set_manager, db_conn, "Ground Truth - Low coverage",
                                   "SV's form genome assembly that are not covered by any alignments", -1).cpp_module.id
    small_caller_run_id = GetCallInserter(parameter_set_manager, db_conn, "Ground Truth - Small",
                                   "SV's form genome assembly that are too small for pacBio reads", -1).cpp_module.id
    one_sided_calls = OneSidedCallsTable(db_conn)

    pack, fm_index, mm_index, query_genomes = load_genomes(query_genome_str, reference_genome)


    cnt_one_sided_jumps = 0
    cnt_low_coverage_jumps = 0
    cnt_remaining_jumps = 0

    idx = 0
    for jumps_, (y_start, query_genome), jump_cov_ in zip(jumps, query_genomes, cov_list):
        query_genome.check()
        for jump, (cov, cov_start, cov_end) in zip(jumps_, jump_cov_):
            jump.id = idx
            assert(jump.num_supp_nt() > 0)
            if cov < min_cov and cov_start >= min_cov and cov_end >= min_cov:
                from_forward = jump.from_forward
                to_forward = jump.to_forward
                if jump.was_mirrored:
                    from_forward = not from_forward
                    to_forward = not to_forward

                from_call = SvCall(jump.from_pos, jump.from_pos, 0, 0, from_forward, from_forward, cov_start,
                                   cov_start * jump.num_supp_nt())
                to_call = SvCall(jump.to_pos, jump.to_pos, 0, 0, to_forward, to_forward, cov_end,
                                 cov_end* jump.num_supp_nt())


                if jump.was_mirrored:
                    to_call.order_id = idx
                    idx += 1
                    from_call.order_id = idx
                    idx += 1
                else:
                    from_call.order_id = idx
                    idx += 1
                    to_call.order_id = idx
                    idx += 1

                if jump.query_from < jump.query_to:
                    if jump.was_mirrored:
                        to_call.inserted_sequence = NucSeq(query_genome, jump.query_from - y_start,
                                                                jump.query_to - y_start)
                        to_call.inserted_sequence.check()
                    else:
                        from_call.inserted_sequence = NucSeq(query_genome, jump.query_from - y_start,
                                                                jump.query_to - y_start)
                        from_call.inserted_sequence.check()

                sv_call_table.insert_call(caller_run_id, from_call) # updates id in from_call
                sv_call_table.insert_call(caller_run_id, to_call) # updates id in to_call

                if jump.was_mirrored:
                    one_sided_calls.insert_calls(to_call, from_call)
                else:
                    one_sided_calls.insert_calls(from_call, to_call)
                cnt_one_sided_jumps += 1
            else:
                if cov < min_cov:
                    cnt_low_coverage_jumps += 1
                else:
                    cnt_remaining_jumps += 1
                call = SvCall(jump.from_pos, jump.to_pos, 0, 0, jump.from_forward, jump.to_forward, cov,
                              cov * jump.num_supp_nt())
                if jump.query_from < jump.query_to:
                    call.inserted_sequence = NucSeq(query_genome, jump.query_from - y_start, jump.query_to - y_start)
                    call.inserted_sequence.check()
                call.order_id = idx
                idx += 1
                call.mirrored = jump.was_mirrored
                if cov < min_cov:
                    sv_call_table.insert_call(low_cov_caller_run_id, call)
                elif max(abs(jump.from_pos - jump.to_pos), abs(jump.query_from - jump.query_to)) < min_size:
                    sv_call_table.insert_call(small_caller_run_id, call)
                else:
                    sv_call_table.insert_call(caller_run_id, call)

    print("Inserted into DB. There were", cnt_one_sided_jumps, "one sided entries,", cnt_low_coverage_jumps,
          "entries with coverage <=", min_cov, "and", cnt_remaining_jumps, "two sided entries with enough coverage" )

    return caller_run_id, None

def render_seeds_2(seeds_1, query_genome, reference_genome, title="seeds"):
    plot = figure(title=title, plot_width=1000, plot_height=1000)

    decorate_plot(plot, query_genome, reference_genome)

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

    plot.xaxis.axis_label = "Reference Genome"
    plot.yaxis.axis_label = "Sequenced Genome"
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

    #make_read_range_table_from_simulated_n_seeds(db_name, seq_id, query_genome, reference_genome)
    out.append(view_coverage(db_name, seq_id, query_genome))

    seeds_n_rects = compute_seeds(query_genome, reference_genome, db_name, seq_id)
    jumps = filter_jumps(seeds_to_jumps(seeds_n_rects, query_genome, reference_genome), reference_genome)
    #jump_coverage_lenient = reads_per_jump(db_name, jumps, seq_id, False, None)
    #jump_coverage = reads_per_jump(db_name, jumps, seq_id, True)
    jump_coverage_strict = reads_per_jump(db_name, jumps, seq_id, True, None)
    #covs = [
    #    (jump_coverage_lenient, "map_q>=0", "red"),
    #    #(jump_coverage, "primary & map_q>=0", "blue"),
    #    (jump_coverage_strict, "primary & map_q>=10", "green")
    #]
    #out.extend(render_reads_per_jump(jumps, covs, query_genome))

    jumps_to_calls_to_db(jumps, jump_coverage_strict, db_name, query_genome, reference_genome)
    #exit()

    out.append(render_seeds_2(seeds_n_rects, query_genome, reference_genome))
    out.append(render_jumps(jumps, jump_coverage_strict, query_genome, reference_genome))
    show(row(out))


# ~/workspace/samtools/samtools view -h minimap2.sorted.bam "chrIX:251682-259682" > minimap2.filtered.bam


# ~/workspace/samtools/samtools view minimap2.filtered.bam | awk '{print length($10), $1}' | sort -n
