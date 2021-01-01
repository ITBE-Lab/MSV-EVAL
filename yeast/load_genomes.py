from MA import *
from MSV import *
from bokeh.plotting import figure, show, reset_output, ColumnDataSource
from bokeh.layouts import column, row
from bokeh.models import FuncTickFormatter
from bokeh.models.tools import HoverTool
from bokeh.io import output_file
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

