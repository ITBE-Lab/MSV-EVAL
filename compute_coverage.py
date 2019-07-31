from sv_jump import *
from MA import *

"""
from intervaltree import Interval, IntervalTree
max_distance = 50

class SvCallPyCoverage:
    def __init__(self, call):
        self.call = call
        # this holds seed intervals for left and right side.
        self.coverage_analysis = {True: [], False: []}
        # ids of contradicting reads
        self.contradicting = []

    def get_from_interval(self):
        return Interval(self.call.from_start, self.call.from_start + self.call.from_size, data=self)

    def get_to_interval(self):
        return Interval(self.call.to_start, self.call.to_start + self.call.to_size, data=self)

    def check_interval(self, start, end, seed, read_id, is_from):
        # check wether this seed matches the from or the to interval (can be both, must be at least one of them)
        if seed.start_ref > end + max_distance or seed.start_ref + seed.size + max_distance < start:
            return False

        # check if the seed belongs to the left or the right region
        is_left = False
        is_right = False
        if seed.start_ref < start:
            is_left = True
        if seed.start_ref + seed.size > end:
            is_right = True

        # we have a contradicting seed
        if is_left and is_right:
            if read_id not in self.contradicting:
                self.contradicting.append(read_id)
        else:
            self.coverage_analysis[is_left].append(seed)

        return True

    def apply_filter(self):
        median_coverage = {True: 0, False: 0}
        for is_left in [True, False]:
            # compute median coverage (via linesweep)
            l = self.coverage_analysis[is_left]
            if len(l) > 0:
                l.sort(key=lambda x: x.start_ref)

                # result list (compressed coverage list [(coverage, num nt the coverage was observed for), ...])
                coverage_list = []
                # seed-start pointer for line sweep
                start_idx = 0
                # seed-end pointer for line sweep
                end_idx = 0
                # current coverage counter
                cov_count = 0
                # last linesweep stop position on reference
                last_pos = l[0].start_ref
                #
                # sweep over seed start and end positions
                # remember the current coverage via counter.
                # for each sweep position: if the sweep position is different than the last one 
                #   (i.e we actually moved forward):
                #       remember the coverage (via counter) and the number of nt the coverage was observed for 
                #       (this pos - last pos)
                # once all positions have been visited sort the compressed list and extract median via while loop
                #
                while start_idx < len(l) or end_idx < len(l):
                    a = l[start_idx].start_ref if start_idx < len(l) else float("inf")
                    b = l[end_idx].start_ref + l[end_idx].size if end_idx < len(l) else float("inf")
                    if a <= b:
                        if a > last_pos:
                            coverage_list.append( (cov_count, a - last_pos) )
                        cov_count += 1
                        last_pos = a
                        start_idx += 1
                    else:
                        if b > last_pos:
                            coverage_list.append( (cov_count, b - last_pos) )
                        cov_count -= 1
                        last_pos = b
                        end_idx += 1
                coverage_list.sort()
                median_cnt = sum(y for x, y in coverage_list) / 2
                idx = 0
                while median_cnt > coverage_list[idx][1]:
                    median_cnt -= coverage_list[idx][1]
                    idx += 1
                #print("cov_list:", coverage_list, idx)
                # we got the median coverage
                median_coverage[is_left] = coverage_list[idx][0]

        # for now just save the median coverage
        self.call.coverage = max(*median_coverage.values(), 1)



    def add_seed(self, seed, read):
        num_found = 0
        if self.check_interval(self.call.from_start, self.call.from_start + self.call.from_size, seed, read.id, True):
            num_found += 1
        if self.check_interval(self.call.to_start, self.call.to_start + self.call.to_size, seed, read.id, False):
            num_found += 1

        # if this did neither overlap the start nor the end interval something is wrong wiht the compute_coverage
        # function below
        assert num_found > 0


def compute_coverage(parameter_set_manager, fm_index, pack, sv_db, caller_id, seq_ids):
    # compute the interval tree of breakpoint-intervals
    calls_from_db = SvCallsFromDb(parameter_set_manager, sv_db, caller_id)
    call_intervals = []
    calls = []
    while calls_from_db.hasNext():
        call = SvCallPyCoverage(calls_from_db.next())
        calls.append(call)
        call_intervals.append(call.get_from_interval())
        call_intervals.append(call.get_to_interval())
    interval_tree = IntervalTree(call_intervals)

    # compute seeds and add them to respective intervals
    for seq_id in seq_ids:
        nuc_seq_getter = AllNucSeqFromSql(parameter_set_manager, sv_db, seq_id)
        seeding_module = BinarySeeding(parameter_set_manager)

        while not nuc_seq_getter.is_finished():
            query = nuc_seq_getter.execute(libMA.ContainerVector())
            segments = seeding_module.execute(fm_index, query)
            seeds = segments.extract_seeds(fm_index, 
                                           parameter_set_manager.by_name("Maximal Ambiguity SV").get(),
                                           parameter_set_manager.by_name("Minimal Seed Size SV").get(),
                                           len(query),
                                           True)

            for seed in seeds:
                for overlap in interval_tree.overlap(seed.start_ref - max_distance,
                                                     seed.start_ref + seed.size + max_distance):
                    overlap.data.add_seed(seed, query)

    for call in calls:
        call.apply_filter()
        sv_db.update_coverage(call.call)
"""

def compute_coverage(parameter_set_manager, fm_index, pack, sv_db, caller_id, seq_ids):
    lock_module = Lock(parameter_set_manager)
    seeding_module = BinarySeeding(parameter_set_manager)
    fm_pledge = Pledge()
    fm_pledge.set(fm_index)
    pack_pledge = Pledge()
    pack_pledge.set(pack)

    for seq_id in seq_ids:
        print("\tconstructing interval tree...")
        coverage_module = libMA.ComputeCoverage(parameter_set_manager, sv_db, caller_id, seq_id)
        print("\tdone constructing interval tree")

        res = VectorPledge()
        # graph for single reads
        for idx in range(parameter_set_manager.get_num_threads()):
            nuc_seq_getter = AllNucSeqFromSql(parameter_set_manager, sv_db, seq_id, idx,
                                              parameter_set_manager.get_num_threads())
            queries_pledge = promise_me(nuc_seq_getter)
            query_pledge = promise_me(lock_module, queries_pledge)
            segments_pledge = promise_me(seeding_module, fm_pledge, query_pledge)
            cov_pledge = promise_me(coverage_module, fm_pledge, queries_pledge, segments_pledge)
            unlock_pledge = promise_me(UnLock(parameter_set_manager, query_pledge), cov_pledge)
            res.append(unlock_pledge)

        # drain all sources
        print("\texecuting graph...")
        res.simultaneous_get( parameter_set_manager.get_num_threads() )
        print("\tdone executing graph")

    print("\tcommitting coverage...")


