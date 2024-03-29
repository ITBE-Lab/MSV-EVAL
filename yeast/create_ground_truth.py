from bokeh.plotting import figure, show, reset_output, ColumnDataSource
from bokeh.layouts import column, row
from bokeh.models import FuncTickFormatter
from bokeh.models.tools import HoverTool
from bokeh.io import output_file
from MA import *
from MSV import *
from sv_util.os_aligners import *
from sv_util.settings import *
from yeast.load_genomes import *
from bisect import bisect_right, bisect_left



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

def jumps_to_calls_to_db(jumps, db_name, query_genome_str, reference_genome, name):
    parameter_set_manager = ParameterSetManager()
    parameter_set_manager.set_selected("SV-PacBio")
    min_size = 200
    db_conn = DbConn({"SCHEMA": {"NAME": db_name}})

    JumpRunTable(db_conn)
    SvCallerRunTable(db_conn)
    sv_call_table = SvCallTable(db_conn)
    #call_per_contig_table = FirstCallPerContigTable(db_conn)
    caller_run_id = GetCallInserter(parameter_set_manager, db_conn, name + " - Large",
                                   "SV's form genome assembly that are too long for Illumina reads", -1).cpp_module.id
    contig_border_caller_run_id = GetCallInserter(parameter_set_manager, db_conn, name + " - Contig Border",
                                   "SV's form genome assembly that are not covered by any alignments", -1).cpp_module.id
    small_caller_run_id = GetCallInserter(parameter_set_manager, db_conn, name + " - Small >=10",
                                   "SV's form genome assembly that are too small for pacBio reads", -1).cpp_module.id
    small_lt_ten_caller_run_id = GetCallInserter(parameter_set_manager, db_conn, name + " - Small <10 ",
                                   "SV's form genome assembly that are too small for pacBio reads", -1).cpp_module.id
    contig_start_caller_run_id = GetCallInserter(parameter_set_manager, db_conn, name + " - Contig Start",
                                   "SV's form genome assembly that are not covered by any alignments", -1).cpp_module.id

    pack, fm_index, mm_index, query_genomes = load_genomes(query_genome_str, reference_genome, parameter_set_manager)


    contig_filter = JumpsFilterContigBorder(parameter_set_manager)

    cnt_small = 0
    cnt_small_lt_10 = 0
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
                # since i do not want to do this in the cpp code i do it here in the py code.
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
                if max(abs(f - t), abs(jump.query_from - jump.query_to)) < 10:
                    sv_call_table.insert_call(small_lt_ten_caller_run_id, call)
                    cnt_small_lt_10 += 1
                else:
                    sv_call_table.insert_call(small_caller_run_id, call)
            else:
                cnt_large += 1
                sv_call_table.insert_call(caller_run_id, call)
            #if is_first:
            #    is_first = False
            #    call_per_contig_table.insert(call.id, query_genome.name, from_forward)

    print("Inserted into DB. There are", cnt_small, "small and", cnt_large, "large entries.", cnt_contig_border,
          "entries are too close to contig borders and are filtered out.", cnt_small_lt_10, "small entries are < 10.")

    return caller_run_id, None

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

    jumps_to_calls_to_db(jumps, db_name, query_genome, reference_genome, "Ground Truth")
    #exit()

    #out.append(render_seeds(seeds_n_rects, query_genome, reference_genome))
    #out.append(render_jumps(jumps, query_genome, reference_genome))
    #show(row(out))

