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
        if dest.score < 10:  # @todo make parameter
            return
        self.destinations.append(dest)

    def get_best_destination(self):
        return max(self.destinations, key=lambda x: x.score)

    def score(self):
        return self.get_best_destination().score


def line_sweep_list(sorted_list, do_functor, get_start, get_end):
    sorted_list.sort(key=lambda x: get_start(x))
    start_it = 0
    end_it = 0

    while start_it < len(sorted_list):
        while end_it < len(sorted_list) and get_end(sorted_list[start_it]) >= get_start(sorted_list[end_it]):
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
                        lambda x: (1 if x.switch_strands else 0, x.ref_pos + fuzziness))

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
    line_sweep_list(sv_jumps, check_sv_jump, lambda x: x.ref_pos, lambda x: x.ref_pos + fuzziness)

    return accepted_sv_jumps


def sv_jumps_to_dict(sv_jumps, accepted_sv_jumps):
    fuzziness = 3
    forw_boxes_data = []
    sw_boxes_data = []
    accepted_boxes_data = []
    for jump in sv_jumps:
        for destination in jump.destinations:
            alpha = 0.08 / math.log(destination.q_distance + 1.5)
            x = [jump.ref_pos - fuzziness, destination.ref_pos -
                 fuzziness, fuzziness*2, fuzziness*2, alpha]
            if destination.switch_strands:
                sw_boxes_data.append(x)
            else:
                forw_boxes_data.append(x)
    for jump in accepted_sv_jumps:
        accepted_boxes_data.append([jump.curr_start, jump.get_best_destination().start, 
                                    jump.curr_end - jump.curr_start,
                                    jump.get_best_destination().end - jump.get_best_destination().start, 0])
    out_dict = {
        "x_offset": 0,
        "panels": [
            {
                "items": [
                    {
                        "type": "box-alpha",
                        "color": "blue",
                        "line_color": None,
                        "line_width": 0,
                        "group": "all_jumps",
                        "data": forw_boxes_data
                    },
                    {
                        "type": "box-alpha",
                        "color": "orange",
                        "line_color": None,
                        "line_width": 0,
                        "group": "all_jumps",
                        "data": sw_boxes_data
                    },
                    {
                        "type": "box-alpha",
                        "color": "#595959",
                        "line_color": "green",
                        "line_width": 3,
                        "group": "accepted_jumps",
                        "data": accepted_boxes_data
                    }
                ],
                "h": 700
            }
        ]
    }
    return out_dict
