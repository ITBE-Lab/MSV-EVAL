from sv_jump import *
import math
from BitVector import BitVector
import sqlite3
from MA import *
from exact_sv_jump_sweep import *
from svCallPy import *
from compute_coverage import *
from analyze_runtimes import AnalyzeRuntimes
import datetime


def sweep_sv_jumps_cpp(parameter_set_manager, sv_db, run_id, ref_size, name, desc, sequencer_ids, pack, fm_index,
                       out_file=None):
    analyze = AnalyzeRuntimes()
    # creates scope so that deconstructor of call inserter is triggered (commits insert transaction)
    def graph():
        print("\tsetting graph up...")
        
        pack_pledge = Pledge()
        pack_pledge.set(pack)

        section_fac = libMA.GenomeSectionFactory(parameter_set_manager, pack)
        lock_module = Lock(parameter_set_manager)
        sweep2 = libMA.ExactCompleteBipartiteSubgraphSweep(parameter_set_manager, sv_db, pack, sequencer_ids[0])

        sv_caller_run_id = sv_db.insert_sv_caller_run(name, desc, run_id)
        filter1 = libMA.FilterLowSupportShortCalls(parameter_set_manager)
        filter2 = libMA.FilterFuzzyCalls(parameter_set_manager)
        filter5 = libMA.FilterDiagonalLineCalls(parameter_set_manager)
        assert len(sequencer_ids) == 1

        res = VectorPledge()
        sections_pledge = promise_me(section_fac) # @note this cannot be in the loop (synchronization!)
        sinks = []
        # graph for single reads
        for _ in range(parameter_set_manager.get_num_threads()):
            # in order to allow multithreading this module needs individual db connections for each thread
            sweep1 = libMA.CompleteBipartiteSubgraphSweep(parameter_set_manager, sv_db, pack, run_id, sequencer_ids[0])
            filter4 = libMA.FilterLowCoverageCalls(parameter_set_manager, sv_db, pack, sequencer_ids[0])
            sink = libMA.BufferedSvCallSink(parameter_set_manager, sv_db, sv_caller_run_id)
            filter3 = libMA.ConnectorPatternFilter(parameter_set_manager, sv_db)
            sinks.append(sink)

            section_pledge = promise_me(lock_module, sections_pledge)
            sweep1_pledge = promise_me(sweep1, section_pledge)
            analyze.register("[0] CompleteBipartiteSubgraphSweep", sweep1_pledge)
            analyze.register("[0.1] CompleteBipartiteSubgraphSweep::init", sweep1, lambda x: x.cpp_module.time_init)
            analyze.register("[0.2] CompleteBipartiteSubgraphSweep::outer_while", sweep1,
                             lambda x: x.cpp_module.time_complete_while - x.cpp_module.time_inner_while)
            analyze.register("[0.3] CompleteBipartiteSubgraphSweep::inner_while", sweep1,
                             lambda x: x.cpp_module.time_inner_while)
            sweep2_pledge = promise_me(sweep2, sweep1_pledge)
            analyze.register("[1] ExactCompleteBipartiteSubgraphSweep", sweep2_pledge)
            #filters

            filter1_pledge = promise_me(filter1, sweep2_pledge)
            analyze.register("[2] FilterLowSupportShortCalls", filter1_pledge)
            filter2_pledge = promise_me(filter2, filter1_pledge)
            analyze.register("[3] FilterFuzzyCalls", filter2_pledge)

            #filter3_pledge = promise_me(filter3, filter2_pledge, pack_pledge) # this filter was off already
            #analyze.register("[4] ConnectorPatternFilter", filter3_pledge)
            #filter3_pledge = promise_me(filter4, filter2_pledge, pack_pledge) # this filter was off already
            #analyze.register("[4] FilterLowCoverageCalls", filter3_pledge)

            filter3_pledge = promise_me(filter5, filter2_pledge)
            analyze.register("[4] FilterDiagonalLineCalls", filter3_pledge)

            write_to_db_pledge = promise_me(sink, filter3_pledge)
            analyze.register("[5] SvCallSink", write_to_db_pledge)
            unlock_pledge = promise_me(UnLock(parameter_set_manager, section_pledge), write_to_db_pledge)
            res.append(unlock_pledge)

        # drain all sources
        print("\texecuting graph...")
        res.simultaneous_get( parameter_set_manager.get_num_threads() )
        print("\tcommiting calls...")
        for sink in sinks:
            sink.cpp_module.commit()
        print("\tdone")

        return sv_caller_run_id

    sv_caller_run_id = graph()
    print("done sweeping")
    print("num calls:", sv_db.get_num_calls(sv_caller_run_id, 0))

    print("overlapping...")
    start = datetime.datetime.now()
    num_combined = libMA.combine_overlapping_calls(parameter_set_manager, sv_db, sv_caller_run_id)
    end = datetime.datetime.now()
    delta = end - start
    analyze.register("[6] combine_overlapping_calls", delta.total_seconds(), lambda x: x)
    print("done overlapping; combined", num_combined, "calls")

    print("computing coverage...")
    start = datetime.datetime.now()
    compute_coverage(parameter_set_manager, fm_index, pack, sv_db, sv_caller_run_id, sequencer_ids)
    end = datetime.datetime.now()
    delta = end - start
    analyze.register("[7] compute_coverage", delta.total_seconds(), lambda x: x)
    print("done computing coverage")

    print("computing coverage index...")
    start = datetime.datetime.now()
    sv_db.add_score_index(sv_caller_run_id)
    end = datetime.datetime.now()
    delta = end - start
    analyze.register("[8] compute_coverage_index", delta.total_seconds(), lambda x: x)
    print("done computing coverage index")

    analyze.analyze(out_file)
    if not out_file is None:
        out_file.write("run_id is " + str(sv_caller_run_id) + "\n")


class YRangeTree:
    def __init__(self, size):
        self.bit_vec = BitVector(size=size)

    def set_to(self, start, end, value):
        assert end > start
        #print("set", start, end, "to", value)
        for x in range(start, end):
            self.bit_vec[x] = value

    def all(self, start, end):
        #print("all", start, end)
        return self.bit_vec[start:end].count_bits() == end - start

    def none(self, start, end):
        return self.bit_vec[start:end].count_bits() == 0

    def find_zero(self, pos):
        assert pos >= 0
        while pos > 0:
            if self.bit_vec[pos] == 0:
                return pos
            else:
                pos -= 1
        return 0

    def find_one(self, pos):
        assert pos >= 0
        while pos > 0:
            if self.bit_vec[pos] == 1:
                return pos
            else:
                pos -= 1
        return 0

    def get_one_intervals_upwards(self, start, end):
        assert end > start
        assert start < len(self.bit_vec)
        ret = []

        if self.none( max(start-1, 0), min(end+1, len(self.bit_vec))):
            return ret

        end = self.find_one(end)+1
        intermediate = self.find_zero(end-1)

        while intermediate >= start:  # while intermediate lower than start
            ret.append((intermediate, end))  # append intermediate interval
            if self.none( max(start-1, 0), intermediate):  # if the rest of the interval is empty
                return ret
            # set end for next interval search
            end = self.find_one(intermediate-1)+1
            intermediate = self.find_zero(end-1)
        ret.append((intermediate, end))  # append the last interval

        return ret


class WarpedBitVector:
    def __init__(self, size, warp_factor=5000, center_strip_up=5000, center_strip_down=1000):
        self.physical_size = math.ceil(   max(0, size - center_strip_up) / warp_factor 
                                        + max(0, size - center_strip_down) / warp_factor)
        self.physical_size += center_strip_up + center_strip_down + 1
        print("physical size:", self.physical_size / (8*10**9), "gb")
        self.bit_vec = YRangeTree(self.physical_size)
        self.warp_factor = warp_factor
        self.center_strip_up = center_strip_up
        self.center_strip_down = center_strip_down
        self.size = size
        assert self.to_new_coord_system_f(self.size - 1, 0) < self.physical_size
        assert self.to_new_coord_system_c(0, self.size - 1) >= 0

    def to_new_coord_system(self, y, x):
        #print("to_new_coord_system", y, x, self.physical_size/2)
        y -= x
        if y > self.center_strip_up:
            y = (y - self.center_strip_up) / self.warp_factor + self.center_strip_up
        if y < -self.center_strip_down:
            y = (y + self.center_strip_down) / self.warp_factor - self.center_strip_down
        # increment by the logical center
        y += self.center_strip_down + int( (self.size - self.center_strip_down) / self.warp_factor)
        #print(y)
        return y

    def to_new_coord_system_f(self, y, x):
        return int(math.floor(self.to_new_coord_system(y, x)))

    def to_new_coord_system_c(self, y, x):
        return int(math.ceil(self.to_new_coord_system(y, x)))

    # returns first coordinate of set interval
    def set_to(self, start, end, value, x_value):
        self.bit_vec.set_to(self.to_new_coord_system_f(start, x_value),
                            self.to_new_coord_system_c(end, x_value), value)
        return self.to_new_coord_system_f(start, x_value)

    def clear_downwards(self, pos):
        while self.bit_vec.bit_vec[pos] == 1:
            self.bit_vec.bit_vec[pos] = 0
            pos += 1

    def all(self, start, end, x_value):
        return self.bit_vec.all(self.to_new_coord_system_f(start, x_value), self.to_new_coord_system_c(end, x_value))

    def none(self, start, end, x_value):
        return self.bit_vec.none(self.to_new_coord_system_f(start, x_value), self.to_new_coord_system_c(end, x_value))

    def find_zero(self, pos, x_value):
        return self.bit_vec.find_zero(self.to_new_coord_system_f(pos, x_value))

    def find_one(self, pos, x_value):
        return self.bit_vec.find_one(self.to_new_coord_system_f(pos, x_value))

    def get_one_intervals_upwards(self, start, end, x_value):
        #print("get_one_intervals_upwards", self.to_new_coord_system_f(
        #    start, x_value), self.to_new_coord_system_c(end, x_value))
        assert end > start
        assert start < self.size
        return self.bit_vec.get_one_intervals_upwards(self.to_new_coord_system_f(start, x_value),
                                                      self.to_new_coord_system_c(end, x_value))


def sweep_sv_jumps(parameter_set_manager, sv_db, run_id, ref_size, name, desc, sequencer_ids, pack, fm_index):
    sweeper = SortedSvJumpFromSql(parameter_set_manager, sv_db, run_id)
    call_inserter = SvCallInserter(sv_db, name, desc, run_id)
    print("creating sweep list...")

    estimated_coverage_list = [0]*len(pack.contigLengths())
    for sequencer_id in sequencer_ids:
        for idx, cnt in enumerate(sv_db.get_num_nts(sequencer_id)):
            estimated_coverage_list[idx] += cnt
    for idx, cnt in enumerate(pack.contigLengths()):
        estimated_coverage_list[idx] /= cnt

    print(name, "- estimated_coverage per contig:", estimated_coverage_list)

    y_range_tree = WarpedBitVector(ref_size)
    cluster_dict = {}
    print("sweeping...")
    #for switch_strand, from_pos, to_start, to_end, is_end, score, jmp in line_sweep_list:
    def sweep_sv_start(sv_jmp):
        #if sv_jmp.to_start() == 0:
        #    print("start:", sv_jmp.to_start(), sv_jmp.to_end(), sv_jmp.from_start(), sv_jmp.from_end())
        #if not sv_jmp.switch_strand_known():
        #    return # @todo
        cluster = SvCallPy(sv_jmp)
        cluster_keys = [x for x, _ in y_range_tree.get_one_intervals_upwards(sv_jmp.to_start(), sv_jmp.to_end(),
                                                                             sv_jmp.from_start_same_strand())]
        #print("cluster_keys", cluster_keys)
        for key in cluster_keys:
            if key not in cluster_dict:
                print("CRITICAL:", key, "not in dict")
                print("bitvec:", y_range_tree.bit_vec.bit_vec[max(0,key-10):min(key+10,len(y_range_tree.bit_vec.bit_vec))])
                assert False
            cluster.join(cluster_dict[key])
            #print("del", key, cluster_dict[key])
            del cluster_dict[key]
        new_key = max(y_range_tree.set_to(sv_jmp.to_start(), sv_jmp.to_end(), 1,
                      sv_jmp.from_start_same_strand()) - 1, 0)
        if len(cluster_keys) > 0:
            new_key = min(new_key, min(cluster_keys))
        #print("insert", new_key, cluster)
        cluster_dict[new_key] = cluster
        if y_range_tree.bit_vec.find_zero(
                y_range_tree.to_new_coord_system_c(sv_jmp.to_end(), sv_jmp.from_start_same_strand()) - 1) != new_key:
            print("inserted and found key do not match",
                  y_range_tree.bit_vec.find_zero(
                      y_range_tree.to_new_coord_system_c(sv_jmp.to_end(), sv_jmp.from_start_same_strand() ) - 1),
                  new_key)
            assert False
    def sweep_sv_end(sv_jmp):
        #if sv_jmp.to_start() == 0:
        #    print("end:", sv_jmp.to_start(), sv_jmp.to_end(), sv_jmp.from_start(), sv_jmp.from_end())
        #if not sv_jmp.switch_strand_known():
        #    return # @todo
        key = y_range_tree.find_zero(sv_jmp.to_start(), sv_jmp.from_start_same_strand())
        if key not in cluster_dict:
            print("CRITICAL:", key, "not in dict")
            print("searched from", y_range_tree.to_new_coord_system_c(
                sv_jmp.to_start(), sv_jmp.from_start_same_strand()))
            print(
                "bitvec:", y_range_tree.bit_vec.bit_vec[max(key-10,0):key+10])
            assert False
        cluster_dict[key].count -= 1
        if len(cluster_dict[key]) <= 0:
            estimated_coverage = min(estimated_coverage_list[pack.seq_id_for_pos(cluster_dict[key].call.from_start)],
                                     estimated_coverage_list[pack.seq_id_for_pos(cluster_dict[key].call.to_start)])
            # check for acceptance:
            # @note these parameters are hardcoded in two locations @todo
            if len(cluster_dict[key].call.supporing_jump_ids) >= max(estimated_coverage/8, 2):
                for accepted_cluster in sweep_sv_call(cluster_dict[key], estimated_coverage):
                    #print("accepting", str(accepted_cluster))
                    call_inserter.insert_call(accepted_cluster.call)
            y_range_tree.clear_downwards(key + 1)
            # @todo are the next two lines correct or not?
            #if key == 0: # fix bug where first bit of cluster is not cleared (not sure wether this is actually correct)
            #    y_range_tree.bit_vec.bit_vec[0] = 0
            #print("del", key, cluster_dict[key])
            del cluster_dict[key]
        #if sv_jmp.to_start() == 0:
        #    print("bitvec:", y_range_tree.bit_vec.bit_vec[0:10])

    while sweeper.has_next_start() and sweeper.has_next_end():
        if sweeper.next_start_is_smaller():
            sweep_sv_start(sweeper.get_next_start())
        else:
            sweep_sv_end(sweeper.get_next_end())
    while sweeper.has_next_start():
        print("something is wierd... (this line should never be reached)")
        sweep_sv_start(sweeper.get_next_start())
    while sweeper.has_next_end():
        sweep_sv_end(sweeper.get_next_end())

    print("done sweeping")
    sv_caller_run_id = call_inserter.sv_caller_run_id
    del call_inserter # trigger deconstructor for call inserter (commits insert transaction)
    print("num calls:", sv_db.get_num_calls(sv_caller_run_id, 0))
    print("filtering low support short calls...")
    num_removed = sv_db.filter_short_edges_with_low_support(sv_caller_run_id, 500, 50)
    print("done filtering; removed", num_removed, "calls")
    print("filtering fuzzy calls calls...")
    num_removed = sv_db.filter_fuzzy_calls(sv_caller_run_id,
                                           parameter_set_manager.by_name("Max Fuzziness Filter").get())
    print("done filtering; removed", num_removed, "calls")
    print("overlapping...")
    num_combined = libMA.combine_overlapping_calls(parameter_set_manager, sv_db, sv_caller_run_id)
    print("done overlapping; combined", num_combined, "calls")
    print("num calls remaining:", sv_db.get_num_calls(sv_caller_run_id, 0))
    print("computing coverage...")
    compute_coverage(parameter_set_manager, fm_index, pack, sv_db, sv_caller_run_id, sequencer_ids)
    print("done computing coverage")


def sv_jumps_to_dict(sv_db, run_ids=None, x=None, y=None, w=None, h=None, only_supporting_jumps=False, min_score=0,
                     max_render = 2000, pack=None):
    forw_boxes_data = []
    unknown_boxes_data_a = []
    unknown_boxes_data_b = []
    sw_boxes_data = []
    accepted_lines_data = []
    plus_data = []
    patch_data = []
    contig_borders = []

    min_ = float("inf")
    max_ = 0

    if run_ids is None:
        run_ids = sv_db.newest_unique_runs( 3 )
    for cnt, run_id in enumerate(run_ids):
        assert sv_db.run_exists(run_id)

        name = sv_db.get_run_name(run_id)
        params = ParameterSetManager()
        if "illumina" in name:
            params.set_selected("SV-Illumina")
        if "pacBio" in name:
            params.set_selected("SV-PacBio")
        if "nanopore" in name:
            params.set_selected("SV-ONT")

        def render_jump(jump):
            nonlocal min_
            nonlocal max_
            xs = [jump.from_start_same_strand() - 0.5,
                    jump.to_start() - 0.5,
                    jump.from_size() + 1,
                    jump.to_size() + 1,
                    jump.num_supp_nt()/1000,
                    "SuppNt: " + str(jump.num_supp_nt()) + " read_id: " + str(jump.read_id) + " q_len: " + 
                    str(jump.query_distance())]
            if jump.switch_strand_known():
                if jump.does_switch_strand():
                    sw_boxes_data.append(xs)
                else:
                    forw_boxes_data.append(xs)
            else:
                if jump.from_known():
                    unknown_boxes_data_a.append(xs)
                else:
                    unknown_boxes_data_b.append(xs)
            f = jump.from_pos
            t = jump.to_pos
            if not jump.from_known():
                f = t
            if not jump.to_known():
                t = f
            min_ = min(f, t, min_)
            max_ = max(f, t, max_)
            if not jump.from_fuzziness_is_rightwards():
                if not jump.to_fuzziness_is_downwards():
                    patch_data.append(
                        [[f - 2.5, f + .5, f + .5],
                        [t - .5, t + 2.5, t - .5]])
                else:
                    patch_data.append(
                        [[f - 2.5, f + .5, f + .5],
                        [t + .5, t - 2.5, t + .5]])
            else:
                if not jump.to_fuzziness_is_downwards():
                    patch_data.append(
                        [[f + 2.5, f - .5, f - .5],
                        [t - .5, t + 2.5, t - .5]])
                else:
                    patch_data.append(
                        [[f + 2.5, f - .5, f - .5],
                        [t + .5, t - 2.5, t + .5]])

        cnt_render = 0
        if only_supporting_jumps:
            calls_from_db = None
            if not None in [x, y, w, h]:
                calls_from_db = SvCallsFromDb(params, sv_db, run_id, x, y, w, h)
            else:
                calls_from_db = SvCallsFromDb(params, sv_db, run_id, min_score)
            while calls_from_db.hasNext():
                cnt_render += 1
                if cnt_render >= max_render:
                    print("hit max_render you wont see the full picture")
                    break
                call = calls_from_db.next()
                if call.num_supp_nt > min_score * call.coverage:
                    for idx in range(len(call.supporing_jump_ids)):
                        render_jump(call.get_jump(idx))
            print("rendered", cnt_render, "calls")
        else:
            sweeper = None
            if not None in [x, y, w, h]:
                sweeper = SortedSvJumpFromSql(params, sv_db, sv_db.get_run_jump_id(run_id), x, y, w, h)
            else:
                sweeper = SortedSvJumpFromSql(params, sv_db, sv_db.get_run_jump_id(run_id))

            cnt_render = 0
            while sweeper.has_next_start():
                cnt_render += 1
                if cnt_render >= max_render:
                    print("hit max_render you wont see the full picture")
                    break
                render_jump(sweeper.get_next_start())
            print("rendered", cnt_render, "jumps")
    if not pack is None:
        for start in pack.contigStarts():
            contig_borders.append([min_, start, max_ - min_, 0])
            contig_borders.append([start, min_, 0, max_ - min_])

    out_dict = {
        "x_offset": 0,
        "panels": [
            {
                "items": [
                    {
                        "type": "line",
                        "color": "lightgray",
                        "group": "contig_borders",
                        "data": contig_borders
                    },
                    {
                        "type": "box-alpha",
                        "color": "blue",
                        "line_color": "blue",
                        "line_width": 3,
                        "group": "all_jumps",
                        "data": forw_boxes_data
                    },
                    {
                        "type": "box-alpha",
                        "color": "grey",
                        "line_color": "grey",
                        "line_width": 3,
                        "group": "all_jumps",
                        "data": unknown_boxes_data_a
                    },
                    {
                        "type": "box-alpha",
                        "color": "yellow",
                        "line_color": "yellow",
                        "line_width": 3,
                        "group": "all_jumps",
                        "data": unknown_boxes_data_b
                    },
                    {
                        "type": "box-alpha",
                        "color": "orange",
                        "line_color": "orange",
                        "line_width": 3,
                        "group": "all_jumps",
                        "data": sw_boxes_data
                    },
                    {
                        "type": "line",
                        "color": "green",
                        "group": "accepted_jumps",
                        "data": accepted_lines_data
                    },
                    {
                        "type": "line",
                        "color": "black",
                        "group": "diagonal",
                        "data": [[min_, min_, max_ - min_, max_ - min_]]
                    },
                    {
                        "type": "patch",
                        "color": "black",
                        "group": "all_jumps",
                        "data": patch_data
                    }
                ],
                "h": 900
            }
        ]
    }

    for cnt, run_id in enumerate(run_ids):
        name = sv_db.get_run_name(run_id)

        params = ParameterSetManager()
        if "illumina" in name:
            params.set_selected("SV-Illumina")
        if "pacBio" in name:
            params.set_selected("SV-PacBio")
        if "nanopore" in name:
            params.set_selected("SV-ONT")

        calls_from_db = None
        if not None in [x, y, w, h]:
            calls_from_db = SvCallsFromDb(params, sv_db, run_id, x, y, w, h)
        else:
            calls_from_db = SvCallsFromDb(params, sv_db, run_id, min_score)
        accepted_boxes_data = []
        accepted_plus_data = []
        cnt_render = 0
        while calls_from_db.hasNext():
            def score(jump):
                if jump.coverage == 0:
                    return ""
                return " score: " + str(jump.num_supp_nt / jump.coverage)
            jump = calls_from_db.next()
            if jump.num_supp_nt > min_score * jump.coverage:
                cnt_render += 1
                if cnt_render >= max_render:
                    print("hit max_render you wont see the full picture")
                    break
                if jump.from_size == 1 and jump.to_size == 1:
                    accepted_plus_data.append([jump.from_start,
                                                jump.to_start,
                                                name + " suppNt: " + str(jump.num_supp_nt) + " cov: " +
                                                str(jump.coverage) + " #reads: " + str(len(jump.supporing_jump_ids)) + 
                                                score(jump)])
                else:
                    accepted_boxes_data.append([jump.from_start - 0.5,
                                                jump.to_start - 0.5,
                                                jump.from_size + 1,
                                                jump.to_size + 1,
                                                0,
                                                name + " suppNt: " + str(jump.num_supp_nt) + " cov: " +
                                                str(jump.coverage) + " #reads: " + str(len(jump.supporing_jump_ids)) + 
                                                score(jump)])
                #if len(jump.l_right) > 0:
                #    accepted_lines_data.append([
                #        jump.right() - 0.5, jump.call.to_start - 0.5,
                #        0, jump.call.to_size + 1
                #    ])
                #if len(jump.l_left) > 0:
                #    accepted_lines_data.append([
                #        jump.left() - 0.5, jump.call.to_start - 0.5,
                #        0, jump.call.to_size + 1
                #    ])
                #if len(jump.l_up) > 0:
                #    accepted_lines_data.append([
                #        jump.call.from_start - 0.5, jump.up() - 0.5,
                #        jump.call.from_size + 1, 0
                #    ])
                #if len(jump.l_down) > 0:
                #    accepted_lines_data.append([
                #        jump.call.from_start - 0.5, jump.down() - 0.5,
                #        jump.call.from_size + 1, 0
                #    ])
        c_list = ["green", "purple", "red", "magenta", "brown", "yellow"]
        sv_call_dict = {
                        "type": "box-alpha",
                        "color": "#595959",
                        "line_color": c_list[cnt],
                        "line_width": 3,
                        "group": "sv_calls",
                        "data": accepted_boxes_data
                    }
        sv_plus_dict = {
                        "type": "plus",
                        "color": c_list[cnt],
                        "group": "sv_calls",
                        "data": accepted_plus_data
                    }
        cnt += 1
        if len(accepted_boxes_data) > 0:
            out_dict["panels"][0]["items"].append(sv_call_dict)
        if len(accepted_plus_data) > 0:
            out_dict["panels"][0]["items"].append(sv_plus_dict)

    #conn = sqlite3.connect(db_name)
    #cur = conn.cursor()
    ## show actual SV crosses
    #cur.execute(
    #    """ SELECT start, end
    #        FROM generated_sv
    #""")
    #for start, end in cur.fetchall():
    #    plus_data.append([start + 0.5, end + 0.5])
    #conn.close()

    return out_dict
