from sv_jump import *
import math
from BitVector import BitVector


class SvJumpCluster:
    def __init__(self, switch_strand, from_pos, to_end, score):
        self.switch_strand = switch_strand
        self.from_start = from_pos
        self.to_end = to_end
        self.count = 1
        self.max_count = 1
        self.score = score
        self.from_end = None
        self.to_start = None

    def join(self, other):
        assert self.switch_strand == other.switch_strand
        self.from_start = min(self.from_start, other.from_start)
        self.to_end = max(self.to_end, other.to_end)
        self.count += other.count
        self.max_count += other.max_count
        self.score += other.score

    def __len__(self):
        return self.count

    def __str__(self):
        return "from " + str(self.from_start) + "~" + str(self.from_end) + " to " + str(self.to_start) + "~" + str(self.to_end) + (" switch_strand" if self.switch_strand else "") + " score " + str(self.score) + " cnt " + str(self.max_count)


class YRangeTree:
    def __init__(self, size):
        self.bit_vec = BitVector(size=size)

    def set_to(self, start, end, value):
        for x in range(start, end):
            self.bit_vec[x] = value

    def all(self, start, end):
        return self.bit_vec[start:end].count_bits() == end - start

    def none(self, start, end):
        return self.bit_vec[start:end].count_bits() == 0

    def find_zero(self, pos):
        jump_dist = 128
        assert pos >= 0
        while jump_dist > 0:
            if self.bit_vec[pos] == 0:
                return pos
            elif jump_dist <= pos and self.all(pos - jump_dist, pos):
                pos -= jump_dist
            else:
                jump_dist = int(jump_dist / 2)
        return pos - 1

    def find_one(self, pos):
        jump_dist = 128
        assert pos >= 0
        while jump_dist > 0:
            if self.bit_vec[pos] == 1:
                return pos
            elif jump_dist <= pos and self.none(pos - jump_dist, pos):
                pos -= jump_dist
            else:
                jump_dist = int(jump_dist / 2)
        return pos - 1

    def get_one_intervals_upwards(self, start, end):
        ret = []

        if self.none(start-1, end+1):
            return ret

        end = self.find_one(end)
        intermediate = self.find_zero(end-1)

        while intermediate > start: # while intermediate lower than start
            ret.append((intermediate, end)) # append intermediate interval
            if self.none(start-1, intermediate-1): # if the rest of the interval is empty
                return ret
            end = self.find_one(intermediate-1) # set end for next interval search
            intermediate = self.find_zero(end-1)
        ret.append((intermediate, end)) # append the last interval

        return ret


def sweep_sv_jumps(parameter_set_manager, conn, sv_jumps, ref_size):
    accepted_sv_jumps = []

    line_sweep_list = []
    print("creating sweep list...")
    for sv_jump in sv_jumps:
        line_sweep_list.append(
            (sv_jump.switch_strands, sv_jump.ref_from_start, sv_jump.ref_to_start, sv_jump.ref_to_end, False, sv_jump.score))
        line_sweep_list.append(
            (sv_jump.switch_strands, sv_jump.ref_from_end, sv_jump.ref_to_start, sv_jump.ref_to_end, True, sv_jump.score))
    print("done creating list")
    print("sweeping...")
    line_sweep_list.sort()
    print("done sorting")

    y_range_tree = YRangeTree(ref_size)
    cluster_dict = {}
    print("sweeping...")
    for switch_strand, from_pos, to_start, to_end, is_end, score in line_sweep_list:
        #print("sweeping", switch_strand, from_pos, to_start, to_end, is_end, score)
        if not is_end:
            cluster = SvJumpCluster(switch_strand, from_pos, to_end, score)
            cluster_keys = [x for x, _ in y_range_tree.get_one_intervals_upwards(to_start, to_end)]
            y_range_tree.set_to(to_start, to_end, 1)
            new_key = to_start - 1
            for key in cluster_keys:
                if key < new_key:
                    new_key = key
                cluster.join(cluster_dict[key])
                del cluster_dict[key]
            #print("insert", new_key, cluster)
            cluster_dict[new_key] = cluster
        else:
            key = y_range_tree.find_zero(to_start)
            cluster_dict[key].count -= 1
            if len(cluster_dict[key]) <= 0:
                # check for acceptance:
                if cluster_dict[key].score >= 0.3 and cluster_dict[key].max_count >= 10:
                    cluster_dict[key].from_end = from_pos
                    cluster_dict[key].to_start = key
                    print("accepting", str(cluster_dict[key]))
                    accepted_sv_jumps.append(cluster_dict[key])
                y_range_tree.set_to(key, cluster_dict[key].to_end, 0)
                #print("delete", key, cluster_dict[key])
                del cluster_dict[key]
    print("done sweeping")


    return accepted_sv_jumps


def sv_jumps_to_dict(sv_jumps, accepted_sv_jumps):
    forw_boxes_data = []
    sw_boxes_data = []
    accepted_boxes_data = []
    for jump in sv_jumps:
        x = [jump.ref_from_start,
             jump.ref_to_start,
             jump.ref_from_end - jump.ref_from_start,
             jump.ref_to_end - jump.ref_to_start,
             jump.score,
             str(jump.read_id)]
        if jump.switch_strands:
            sw_boxes_data.append(x)
        else:
            forw_boxes_data.append(x)
    for jump in accepted_sv_jumps:
        accepted_boxes_data.append([jump.from_start,
                                    jump.to_start,
                                    jump.from_end - jump.from_start,
                                    jump.to_end - jump.to_start,
                                    0,
                                    str(jump.score)])
    out_dict = {
        "x_offset": 0,
        "panels": [
            {
                "items": [
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
                        "color": "orange",
                        "line_color": "orange",
                        "line_width": 3,
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
