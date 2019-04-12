from sv_jump import *
import math
from BitVector import BitVector
import sqlite3


class SvJumpCluster:
    def __init__(self, switch_strand, from_pos, to_start, to_end, score):
        self.switch_strand = switch_strand
        self.from_start = from_pos
        self.to_end = to_end
        self.count = 1
        self.max_count = 1
        self.score = score
        self.from_end = None
        self.to_start = to_start

    def join(self, other):
        if not self.switch_strand == other.switch_strand:
            print("CRITICAL CANNOT JOIN...")
            return
        self.from_start = min(self.from_start, other.from_start)
        self.to_end = max(self.to_end, other.to_end)
        self.to_start = min(self.to_start, other.to_start)
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
    def __init__(self, size, warp_factor=5000, center_strip_up=5000, center_strip_down=100):
        self.physical_size = int(
            (size*2 - (center_strip_up + center_strip_down)) / warp_factor) + (center_strip_up + center_strip_down)
        print("physical size:", self.physical_size / (8*10**9), "gb")
        self.bit_vec = YRangeTree(self.physical_size)
        self.warp_factor = warp_factor
        self.center_strip_up = center_strip_up
        self.center_strip_down = center_strip_down
        self.size = size

    def to_new_coord_system(self, y, x):
        #print("to_new_coord_system", y, x, self.physical_size/2)
        y -= x
        if y > self.center_strip_up:
            y = (y - self.center_strip_up) / self.warp_factor + self.center_strip_up
        if y < -self.center_strip_down:
            y = (y + self.center_strip_down) / self.warp_factor - self.center_strip_down
        y += self.physical_size / 2
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
        return self.bit_vec.get_one_intervals_upwards(self.to_new_coord_system_f(start, x_value),
                                                      self.to_new_coord_system_c(end, x_value))


def sweep_sv_jumps(parameter_set_manager, conn, sv_jumps, ref_size):
    accepted_sv_jumps = []

    line_sweep_list = []
    print("creating sweep list...")
    for sv_jump in sv_jumps:
        assert sv_jump.ref_to_start < sv_jump.ref_to_end
        line_sweep_list.append(
            (sv_jump.switch_strands, sv_jump.ref_from_start, sv_jump.ref_to_start, sv_jump.ref_to_end, False, sv_jump.score, sv_jump.ref_from_start))
        line_sweep_list.append(
            (sv_jump.switch_strands, sv_jump.ref_from_end, sv_jump.ref_to_start, sv_jump.ref_to_end, True, sv_jump.score, sv_jump.ref_from_start))
    print("done creating list")
    print("sweeping...")
    line_sweep_list.sort()
    print("done sorting")

    y_range_tree = WarpedBitVector(ref_size)
    cluster_dict = {}
    print("sweeping...")
    for switch_strand, from_pos, to_start, to_end, is_end, score, from_start in line_sweep_list:
        #print("sweeping", switch_strand, from_pos,
        #      to_start, to_end, is_end, score)
        if not is_end:
            cluster = SvJumpCluster(
                switch_strand, from_pos, to_start, to_end, score)
            cluster_keys = [
                x for x, _ in y_range_tree.get_one_intervals_upwards(to_start, to_end, from_pos)]
            #print("cluster_keys", cluster_keys)
            for key in cluster_keys:
                if key not in cluster_dict:
                    print("CRITICAL:", key, "not in dict")
                    print(
                        "bitvec:", y_range_tree.bit_vec.bit_vec[key-10:key+10])
                    assert False
                cluster.join(cluster_dict[key])
                #print("del", key, cluster_dict[key])
                del cluster_dict[key]
            new_key = max(y_range_tree.set_to(to_start, to_end, 1, from_pos) - 1, 0)
            if len(cluster_keys) > 0:
                new_key = min(new_key, min(cluster_keys))
            #print("insert", new_key, cluster)
            cluster_dict[new_key] = cluster
            if y_range_tree.bit_vec.find_zero(y_range_tree.to_new_coord_system_c(to_end, from_pos) - 1) != new_key:
                print("inserted and found key do not match", 
                      y_range_tree.bit_vec.find_zero(y_range_tree.to_new_coord_system_c(to_end, from_pos) - 1),
                      new_key)
                assert False
        else:
            key = y_range_tree.find_zero(to_start, from_start)
            if key not in cluster_dict:
                print("CRITICAL:", key, "not in dict")
                print("searched from", y_range_tree.to_new_coord_system_c(
                    to_start, from_start))
                print(
                    "bitvec:", y_range_tree.bit_vec.bit_vec[key-10:key+10])
                assert False
            cluster_dict[key].count -= 1
            if len(cluster_dict[key]) <= 0:
                # check for acceptance:
                if cluster_dict[key].score >= 0.3 and cluster_dict[key].max_count >= 10:
                    cluster_dict[key].from_end = from_pos
                    print("accepting", str(cluster_dict[key]))
                    accepted_sv_jumps.append(cluster_dict[key])
                y_range_tree.clear_downwards(key + 1)
                #print("del", key, cluster_dict[key])
                del cluster_dict[key]
    print("done sweeping")

    return accepted_sv_jumps


def sv_jumps_to_dict(sv_jumps, accepted_sv_jumps, db_name):
    forw_boxes_data = []
    sw_boxes_data = []
    accepted_boxes_data = []
    plus_data = []
    patch_data = []
    for jump in sv_jumps:
        x = [jump.ref_from_start - 0.5,
             jump.ref_to_start - 0.5,
             jump.ref_from_end - jump.ref_from_start + 1,
             jump.ref_to_end - jump.ref_to_start + 1,
             jump.score,
             str(jump.read_id)]
        if jump.switch_strands:
            sw_boxes_data.append(x)
        else:
            forw_boxes_data.append(x)
        if jump.fuzziness_from_dir == "left":
            if jump.fuzziness_to_dir == "up":
                patch_data.append(
                    [[jump.x - 2.5, jump.x + .5, jump.x + .5], [jump.y - .5, jump.y + 2.5, jump.y - .5]])
            if jump.fuzziness_to_dir == "down":
                patch_data.append(
                    [[jump.x - 2.5, jump.x + .5, jump.x + .5], [jump.y + .5, jump.y - 2.5, jump.y + .5]])
        else:
            if jump.fuzziness_to_dir == "up":
                patch_data.append(
                    [[jump.x + 2.5, jump.x - .5, jump.x - .5], [jump.y - .5, jump.y + 2.5, jump.y - .5]])
            if jump.fuzziness_to_dir == "down":
                patch_data.append(
                    [[jump.x + 2.5, jump.x - .5, jump.x - .5], [jump.y + .5, jump.y - 2.5, jump.y + .5]])
    for jump in accepted_sv_jumps:
        accepted_boxes_data.append([jump.from_start - 0.5,
                                    jump.to_start - 0.5,
                                    jump.from_end - jump.from_start + 1,
                                    jump.to_end - jump.to_start + 1,
                                    0,
                                    str(jump.score)])
    conn = sqlite3.connect(db_name)
    cur = conn.cursor()
    # show actual SV crosses
    cur.execute(
        """ SELECT start, end
            FROM generated_sv
    """)
    for start, end in cur.fetchall():
        plus_data.append([start + 0.5, end + 0.5])
    conn.close()

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
                    },
                    {
                        "type": "patch",
                        "color": "black",
                        "group": "all_jumps",
                        "data": patch_data
                    },
                    {
                        "type": "plus",
                        "color": "red",
                        "group": "actual_sv",
                        "data": plus_data
                    }
                ],
                "h": 900
            }
        ]
    }
    return out_dict
