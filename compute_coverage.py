from sv_jump import *
from MA import *
from intervaltree import Interval, IntervalTree

max_distance = 25

class SvCallPyCoverage:
    def __init__(self, call):
        self.call = call
        # this holds seed intervals for left and right side.
        self.coverage_analysis = {True: [], False: []}
        self.contradicting = []

    def get_from_interval(self):
        return Interval(self.call.from_start(), self.call.from_size(), self)

    def get_to_interval(self):
        return Interval(self.call.to_start(), self.call.to_size(), self)

    def check_interval(self, start, end, seed, supporting_read, is_from):
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
            self.contradicting.append(seed)
        else:
            self.coverage_analysis[is_left].append(seed)

        return True

    def apply_filter(self):
        self.contradicting.sort(key=lambda x: x.start_ref)
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
                while start_idx < len(l) of end_idx < len(l):
                    a = l[start_idx].start_ref if start_idx < len(l) else float("inf")
                    b = l[end_idx].start_ref + l[end_idx].size if end_idx < len(l) else float("inf")
                    if a <= b:
                        if a > last_pos:
                            coverage_list.append( (cov_count, a - last_pos) )
                        cov_count += 1
                        last_pos = a
                    else:
                        if b > last_pos:
                            coverage_list.append( (cov_count, b - last_pos) )
                        cov_count -= 1
                        last_pos = b
                coverage_list.sort()
                median_cnt = sum(y for x, y in coverage_list) / 2
                idx = 0
                while median_cnt > coverage_list[idx][1]:
                    median_cnt -= coverage_list[idx][1]
                    idx += 1
                # we got the median coverage
                median_coverage = coverage_list[idx][0]

                #@todo continue here tomorrow



    def add_seed(self, seed, read):
        supporting_read = read.id in self.call.supporing_jump_ids:

        num_found = 0
        if check_interval(self.call.from_start(), self.call.from_size(), seed, supporting_read, True):
            num_found += 1
        if check_interval(self.call.to_start(), self.call.to_size(), seed, supporting_read, False):
            num_found += 1

        # if this did neither overlap the start nor the end interval something is wrong wiht the compute_coverage
        # function below
        assert num_found > 0


def compute_coverage(parameter_set_manager, fm_index, pack, sv_db, caller_id, seq_ids):
    # compute the interval tree of breakpoint-intervals
    calls_from_db = SvCallsFromDb(parameter_set_manager, sv_db, caller_id)
    call_intervals = []
    while calls_from_db.has_next():
        call = SvCallPyCoverage(calls_from_db.next())
        call_intervals.append(call.get_from_interval())
        call_intervals.append(call.get_to_interval())
    interval_tree = IntervalTree(call_intervals)

    # compute seeds and add them to respective intervals
    for seq_id in seq_ids:
        nuc_seq_getter = AllNucSeqFromSql(parameter_set_manager, sv_db, seq_id)
        seeding_module = BinarySeeding(parameter_set_manager)

        while not nuc_seq_getter.is_finished():
            query = nuc_seq_getter.execute()
            seeds = seeding_module.execute(fm_index, query)

            for seed in seeds:
                for overlap in interval_tree.overlap(seed.start_ref - max_distance,
                                                     seed.start_ref + seed.size + max_distance):
                    overlap.data.add_seed(seed, query)

    for call_interval in call_intervals:
        call_interval.apply_filter()


