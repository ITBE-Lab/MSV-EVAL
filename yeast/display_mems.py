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


    #print("seq", pack.extract_from_to(185526, 187197))
    #exit()

    return pack, fm_index, mm_index, ret_query_genome

def compute_seeds(seeder, query_genome, reference_genome, ambiguity=10000, seed_filter=None):
    param = ParameterSetManager()
    param.by_name("Max Size Reseed").set(2000)

    pack, fm_index, mm_index, query_genomes = load_genomes(query_genome, reference_genome)
    mm_index.set_max_occ(ambiguity)

    ret = []
    for y_start, query_genome in query_genomes:
        seeds = seeder.execute(mm_index, query_genome, pack)
        mems = SeedLumping(param).execute(seeds, query_genome, pack)
        if not seed_filter is None:
            mems_2 = seed_filter(mems, str(query_genome))
            print("filtered from", len(mems), "down to",len(mems_2))
            mems = mems_2
        soc_filter = SoCFilter(mems)
        filtered_2 = Seeds(soc_filter.pop_all_larger_than())

        jumps_from_seeds = SvJumpsFromSeeds(param, pack)
        reseeding = RecursiveReseeding(param, pack)
        helper_ret = jumps_from_seeds.cpp_module.execute_helper(filtered_2, pack, query_genome)
        #reseeded_mems = helper_ret.seeds
        #layer_of_seeds = helper_ret.layer_of_seeds
        rectangles = helper_ret.rectangles
        parlindromes = helper_ret.parlindrome
        #fill_of_rectangles = helper_ret.rectangles_fill
        #seed_sample_sizes = helper_ret.rectangle_ambiguity
        #rectangle_used_dp = helper_ret.rectangle_used_dp
        ret.append((y_start, reseeding.execute(filtered_2, pack, query_genome), rectangles))
    return ret

def seeds_to_filter_sections(seeds, add_size=50, min_nt=25):
    intervals = []
    for seed in seeds:
        if seed.on_forward_strand:
            intervals.append((seed.start_ref, seed.start_ref + seed.size, seed.size))
        else:
            intervals.append((seed.start_ref - seed.size, seed.start_ref, seed.size))
    intervals.sort()
    max_len_idx = 0
    for idx, interval in enumerate(intervals):
        if interval[2] > intervals[max_len_idx][2]:
            max_len_idx = idx
    del intervals[max_len_idx]

    merged_intervals = [[*intervals[0]]]
    for start, end, nt in intervals[1:]:
        if start <= merged_intervals[-1][1] + add_size:
            merged_intervals[-1][1] = end
            merged_intervals[-1][2] += nt
        else:
            merged_intervals.append([start, end, nt])

    ret_intervals = [(start, end) for start,end,num in merged_intervals if num > min_nt]

    return ret_intervals

def str_to_k_mer(s, k=18):
    for idx in range(1 + len(s) - k):
        sec = s[idx:idx+k]
        yield sec
        nuc_seq = NucSeq(sec)
        nuc_seq.reverse()
        nuc_seq.complement()
        yield str(nuc_seq)

def filter_k_mer_set(genome, max_cnt=1):
    k_mers = {}
    for k_mer in str_to_k_mer(genome):
        if not k_mer in k_mers:
            k_mers[k_mer] = 0 
        k_mers[k_mer] += 1

    if False:
        with open("tmp_k_mers_assembly.txt", "w") as out_file:
            for k_mer, cnt in k_mers.items():
                if cnt > 1:
                    out_file.write(k_mer)
                    out_file.write("\t")
                    out_file.write(str(cnt))
                    out_file.write("\n")
        exit()

    return set([k_mer for k_mer, cnt in k_mers.items() if cnt > max_cnt])

def load_filter_set(file_name="tmp_k_mers.txt", t=300):
    s = set()
    with open(file_name, "r") as in_file:
        for line in in_file.readlines():
            seq, cnt = line.split("\t")
            if int(cnt) > t:
                s.add(seq)
    print("query filter k-mers:", len(s))
    return s

k_mers_filter = load_filter_set()
def filter_by_k_mer_set(seeds, query):
    #k_mers = filter_k_mer_set(query)
    ret = []
    for seed in seeds:
        q_str = query[seed.start:seed.start + seed.size]
        if 'n' in q_str or 'N' in q_str:
            continue
        f = False
        #print("q", q_str)
        for k_mer in str_to_k_mer(q_str):
            #print(k_mer)
            if k_mer in k_mers_filter:
                f = True
                break
        if not f:
            ret.append(seed)
    return Seeds(ret)

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

class SoCFilter:
    def __init__(self, seeds):
        self.seeds = [s for s in seeds]
        self.max_delta_dist = 50
        self.socs = [] # start_idx, end_idx, num_nt
        self.__soc_sweep()
        self.__seperate_socs()

    def delta(self, seed):
        if seed.on_forward_strand:
            return seed.start_ref - seed.start
        return seed.start_ref + seed.start + seed.size

    def delta_dist(self, seed_a, seed_b):
        if seed_a.on_forward_strand != seed_b.on_forward_strand:
            return False
        return abs(self.delta(seed_a) - self.delta(seed_b)) < self.max_delta_dist

    def __soc_sweep(self):
        self.seeds.sort(key=lambda x:(x.on_forward_strand, x.start_ref-x.start))
        start_idx = 0
        end_idx = 0
        num_nt = 0
        while start_idx < len(self.seeds):
            while end_idx < len(self.seeds) and self.delta_dist(self.seeds[start_idx], self.seeds[end_idx]):
                num_nt += self.seeds[end_idx].size
                end_idx += 1
            if len(self.socs) > 0 and self.delta_dist(self.seeds[start_idx], self.seeds[self.socs[-1][1]-1]):
                if len(self.socs) > 0 and num_nt > self.socs[-1][2]:
                    del self.socs[-1]
                else:
                    num_nt -= self.seeds[start_idx].size
                    start_idx += 1
                    continue
            self.socs.append((start_idx, end_idx, num_nt))
            num_nt -= self.seeds[start_idx].size
            start_idx += 1

    def __seperate_socs(self):
        new_socs = []
        for start_idx, end_idx, _ in self.socs:
            seeds = self.seeds[start_idx:end_idx]
            seeds.sort(key=lambda x:x.start)
            max_q = -self.max_delta_dist - 10
            for seed in seeds:
                if seed.start >= max_q + self.max_delta_dist:
                    new_socs.append([])
                new_socs[-1].append(seed)
                max_q = max(max_q, seed.start + seed.size)
        new_socs.sort(key=lambda x: sum([s.size for s in x]))
        self.socs = new_socs

    def pop(self):
        ret = self.socs[-1]
        del self.socs[-1]
        return ret

    def pop_all_larger_than(self, num_nt=150):
        ret = []
        while len(self.socs) > 0:
            p = self.pop()
            if sum([s.size for s in p]) < num_nt:
                return ret
            for s in p:
                ret.append(s)
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
        ret.append(filtered_jumps)
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
            xs.append(f)
            ys.append(t)
            ms.append(jump.was_mirrored)
            cs.append(colors[jump.from_forward][jump.to_forward])
            ids.append(idx)
            idx+=1
    plot = figure(title="entries", plot_width=1000, plot_height=1000)
    decorate_plot(plot, reference_genome, reference_genome, True)
    plot.x(x="xs", y="ys", color="cs", line_width=4, source=ColumnDataSource(data={
                    "xs":xs, "ys":ys, "cs":cs, "ms":ms, "ids": ids
                }))
    plot.add_tools(HoverTool(tooltips=[("x, y", "@xs, @ys"),("id", "@ids"),("color, mirrored", "@cs, @ms")]))
    return plot

def jumps_to_calls_to_db(jumps, db_name, query_genome_str, reference_genome):
    db_conn = DbConn({"SCHEMA": {"NAME": db_name}})

    JumpRunTable(db_conn)
    SvCallerRunTable(db_conn)
    get_inserter = GetCallInserter(ParameterSetManager(), db_conn, "Ground Truth",
                                   "SV's form genome assembly", -1)
    pool = PoolContainer(1, db_name)
    sv_inserter = get_inserter.execute(pool)
    
    pack, fm_index, mm_index, query_genomes = load_genomes(query_genome_str, reference_genome)

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
    param = ParameterSetManager()

    seeds_n_rects = compute_seeds(MinimizerSeeding(param), query_genome, reference_genome, 2, filter_by_k_mer_set)
    jumps = seeds_to_jumps(seeds_n_rects, query_genome, reference_genome)
    jumps_to_calls_to_db(jumps, "UFRJ50816", query_genome, reference_genome) #
    #exit()

    aligner_seeds = None
    #aligner_seeds = run_aligner(query_genome, reference_genome)
    #aligner_seeds = [(x,y,[]) for x,y in aligner_seeds]

    out = []
    out.append(render_seeds_2(seeds_n_rects, aligner_seeds, query_genome, reference_genome))
    out.append(render_jumps(jumps, query_genome, reference_genome))
    show(row(out))