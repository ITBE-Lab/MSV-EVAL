from sv_jump import *
import math


class AcceptedSvJump:
    def __init__(self, curr_jump, read_id_set):
        self.curr_start = min(x.ref_pos for x in curr_jump)
        self.curr_end = max(x.ref_pos for x in curr_jump)
        self.destinations = []

    def overlaps(self, other):
        if self.curr_end >= other.curr_start and self.curr_start <= other.curr_end:
            return True
        return False

    def add_destination(self, curr_destinations, dest_read_id_set):
        class Destination:
            def __init__(self, curr_destinations, dest_read_id_set):
                self.start = min(x.ref_pos for x in curr_destinations)
                self.end = max(x.ref_pos for x in curr_destinations)

                self.dist_list = [x.q_distance for x in curr_destinations]
                self.case_list = [x.case for x in curr_destinations]
                self.score = self.score_for_distance_list(self.dist_list)
                self.supporting = len(dest_read_id_set)
                self.switch_strands = curr_destinations[0].switch_strands

            def score_for_distance_list(self, dist_list):
                score = 0
                for x in dist_list:
                    score += 1 / math.log(x + 1.5)
                return score

        dest = Destination(curr_destinations, dest_read_id_set)
        if dest.end + 3 >= self.curr_start and dest.start <= self.curr_end + 3:
            return
        for x in self.destinations:
            if x.start <= dest.end and x.end >= dest.start:
                return
        if dest.score < 10: # @todo make parameter
            return
        self.destinations.append(dest)

    def get_best_destination(self):
        return max(self.destinations, key=lambda x: x.score)

    def score(self):
        return self.get_best_destination().score


def line_sweep_list(sorted_list, do_functor, key, fuzziness):
    # returns whether a and b are closer than d
    # works on tuples...
    def is_closer(a, b, d):
        if isinstance(a, tuple) and isinstance(a, tuple) and isinstance(d, tuple):
            assert len(a) == len(b) == len(d)
            for a_i, b_i, d_i in zip(a,b,d):
                if b_i - d_i < a_i:
                    return True
                elif b_i - d_i > a_i:
                    return False
            return True
        return b - d <= a
    sorted_list.sort(key=lambda x: key(x))
    start_it = 0
    end_it = 0

    while start_it < len(sorted_list):
        while end_it < len(sorted_list) and is_closer(key(sorted_list[start_it]), key(sorted_list[end_it]), fuzziness):
            end_it += 1
        do_functor(sorted_list[start_it:end_it])
        start_it += 1


def sweep_sv_jumps(parameter_set_manager, conn, sv_jumps):
    accepted_sv_jumps = []
    fuzziness = parameter_set_manager.get_selected().by_name("fuzziness").get()
    min_coverage = parameter_set_manager.get_selected().by_name("min coverage").get()
    cursor = conn.cursor()
    cursor.execute("DROP TABLE IF EXISTS sv_line")
    cursor.execute(
        "CREATE TABLE sv_line (id INTEGER PRIMARY KEY, start INTEGER, end INTEGER)")
    cursor.execute("DROP TABLE IF EXISTS sv_jump")
    cursor.execute("""
        CREATE TABLE sv_jump
           (id PRIMARY KEY,
            sv_line_id INTEGER,
            start INTEGER,
            end INTEGER,
            switch_strand BOOL)
        """)
    last_line_id = None

    def check_sv_jump(curr_jumps):
        nonlocal last_line_id
        read_id_set = set()
        for curr_jump in curr_jumps:
            read_id_set.add(curr_jump.read_id)
        if len(read_id_set) < min_coverage:
            return
        # extract all destinations
        dest_list = []
        for curr_jump in curr_jumps:
            dest_list.extend(curr_jump.destinations)

        curr = AcceptedSvJump(curr_jumps, read_id_set)

        # @note this is a temp filter for testing...
        if curr.curr_start > 8000000:
            return

        # and sweep again
        def accept_functor(curr_destinations):
            dest_read_id_set = set()
            for curr_dest in curr_destinations:
                dest_read_id_set.add(curr_dest.parent.read_id)
            if len(dest_read_id_set) < min_coverage:
                return

            curr.add_destination(curr_destinations, dest_read_id_set)
        line_sweep_list(dest_list,
                        accept_functor,
                        lambda x: (1 if x.switch_strands else 0, x.ref_pos),
                        (0, fuzziness))

        if len(curr.destinations) == 0:
            return
        # if this and the last jump overlap keep the better scored one merely
        if len(accepted_sv_jumps) > 0 and accepted_sv_jumps[-1].overlaps(curr):
            if curr.score() <= accepted_sv_jumps[-1].score():
                return
            else:
                print("!!! replacing last sv jump !!!")
                cursor.execute(
                    "DELETE FROM sv_jump WHERE sv_line_id == ?", (last_line_id,))
                cursor.execute(
                    "DELETE FROM sv_line WHERE id == ?", (last_line_id,))
                del accepted_sv_jumps[-1]
        accepted_sv_jumps.append(curr)
        cursor.execute("INSERT INTO sv_line VALUES (NULL, ?, ?)",
                       (curr.curr_start, curr.curr_end))
        last_line_id = cursor.lastrowid
        dest = curr.get_best_destination()
        cursor.execute("INSERT INTO sv_jump VALUES (NULL, ?, ?, ?, ?)",
                       (last_line_id, dest.start, dest.end, dest.switch_strands))

        print("\njump from:", curr.curr_start)
        print("", "to", "score", "sw_st", "dist_list", sep="\t")
        for x in sorted(curr.destinations, key=lambda x: x.score, reverse=True):
            print("", x.start, x.score, x.switch_strands,
                  x.dist_list, x.case_list, sep="\t")

    # line sweep
    line_sweep_list(sv_jumps, check_sv_jump, lambda x: x.ref_pos, fuzziness)

    return accepted_sv_jumps
