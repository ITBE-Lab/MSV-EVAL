from bokeh.plotting import figure, show, reset_output, ColumnDataSource
from MA import *
from sv_util.os_aligners import *
from bisect import bisect_right

global_prefix = "/MAdata/"
genome_dir = global_prefix + "genome/yeasts/"

def load_genomes(query_genome, reference_genome):
    pack = Pack()
    pack.load(genome_dir + reference_genome + "/ma/genome")
    fm_index = FMIndex()
    fm_index.load(genome_dir + reference_genome + "/ma/genome")
    mm_index = MinimizerIndex(ParameterSetManager(), libMA.util.StringVector([str(pack.extract_forward_strand())]), libMA.util.StringVector(["chrVII"]))
    file_reader = FileReader(ParameterSetManager())
    f_stream = FileStreamFromPath(genome_dir + query_genome + "/fasta/genome.fna")
    query_genome = NucSeq()
    while not f_stream.eof():
        query_genome.append(str(file_reader.execute(f_stream)))
    query_genome.name = "sequence0"

    return pack, fm_index, mm_index, query_genome

def compute_seeds(seeder, query_genome, reference_genome, ambiguity=10000):
    pack, fm_index, mm_index, query_genome = load_genomes(query_genome, reference_genome)
    mm_index.set_max_occ(ambiguity)

    print("a")
    seeds = seeder.execute(mm_index, query_genome, pack)
    print("b")
    mems = SeedLumping(ParameterSetManager()).execute(seeds, query_genome, pack)
    print("c")
    filtered_mems = MinLength(ParameterSetManager(), 20).execute(mems)
    print("d")
    #seeds = segments.extract_seeds(fm_index, 10000, 18, len(query_genome), True)
    return filtered_mems

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

def str_to_k_mer(s, k=16):
    for idx in enumerate(len(s) - k):
        yield s[idx:idx+k]

def filter_k_mer_set(genome, k=16, max_cnt=5):
    k_mers = {}
    pack = Pack()
    pack.load(genome_dir + genome + "/ma/genome")
    gen_str = pack.extract_forward_strand()
    for k_mer in str_to_k_mer(gen_str):
        if not k_mer in k_mers:
            k_mers[k_mer] = 0
        k_mers[k_mer] += 1

    return set([k_mer for k_mer, cnt in k_mers.items() if cnt < max_cnt])

def filter_by_k_mer_set(seeds, query, ref, k_mers):
    ret = []
    for seed in seeds:
        q_str = query[seed.start:seed.start + seed.size]
        if seed.on_forward_strand:
            r_str = ref[seed.start_ref:seed.start_ref + seed.size]
        else:
            r_str = ref[seed.start_ref - seed.size :seed.start_ref]
        f = False
        for k_mer in str_to_k_mer(q_str):
            if k_mer in k_mers:
                f = True
        for k_mer in str_to_k_mer(r_str):
            if k_mer in k_mers:
                f = True
        if not f:
            ret.append(seed)

def run_aligner(query_genome_str, reference_genome):
    pack, fm_index, mm_index, query_genome = load_genomes(query_genome_str, reference_genome)

    with open("sequenced.fasta", "w") as out_file:
        out_file.write(">sequence0\n")
        out_file.write(str(query_genome))
        out_file.write("\n")


    read_set = {"technology":"pb", "name":"test", "fasta_file":"sequenced.fasta"}
    json_dict = {"reference_path":genome_dir + reference_genome}
    sam_file_path = "mm2.sam"
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
    return seeds_list

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

    def pop_all_larger_than(self, num_nt=200):
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
    show(plot)

#query_genome = "knowlesiStrain"
#query_genome = "UFRJ50816-chrVII-section"
query_genome = "UFRJ50816"
#query_genome = "YPS138-chrVII-section"
#reference_genome = "YPS138-chrVII-section"
reference_genome = "YPS138"
#reference_genome = "vivax"
if __name__ == "__main__":
    param = ParameterSetManager()
    #seeds = compute_seeds(MinimizerSeeding(param), reference_genome, reference_genome)
    #intervals = seeds_to_filter_sections(seeds)
    #print("#intervals:", len(intervals))
    seeds_y = compute_seeds(MinimizerSeeding(param), query_genome, query_genome)
    intervals_y = seeds_to_filter_sections(seeds_y)
    seeds_2 = compute_seeds(MinimizerSeeding(param), query_genome, reference_genome, 2)
    seed_filter = SeedsFilter(intervals_2=intervals_y)
    print("e")
    filtered_1 = [s for s in seeds_2 if not seed_filter.seed_filtered(s)]
    print("f")
    soc_filter = SoCFilter(filtered_1)
    print("g")
    filtered_2 = soc_filter.pop_all_larger_than()
    print("h")
    render_seeds(filtered_2)
    #render_seeds(run_aligner(query_genome, reference_genome))