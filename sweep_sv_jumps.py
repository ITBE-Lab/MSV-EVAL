from sv_jump import *
import math
from BitVector import BitVector
import sqlite3
from MA import *
from exact_sv_jump_sweep import *
from svCallPy import *

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
        self.physical_size = int((size - center_strip_up) / warp_factor + (size - center_strip_down) / warp_factor)
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


def sweep_sv_jumps(parameter_set_manager, sv_db, run_id, ref_size):
    sweeper = SortedSvJumpFromSql(sv_db, run_id)
    call_inserter = SvCallInserter(sv_db, run_id)
    print("creating sweep list...")

    y_range_tree = WarpedBitVector(ref_size)
    cluster_dict = {}
    print("sweeping...")
    #for switch_strand, from_pos, to_start, to_end, is_end, score, jmp in line_sweep_list:
    def sweep_sv_start(sv_jmp):
        cluster = SvCallPy(sv_jmp)
        cluster_keys = [x for x, _ in y_range_tree.get_one_intervals_upwards(sv_jmp.to_start(), sv_jmp.to_end(),
                                                                             sv_jmp.from_start_same_strand())]
        #print("cluster_keys", cluster_keys)
        for key in cluster_keys:
            if key not in cluster_dict:
                print("CRITICAL:", key, "not in dict")
                print("bitvec:", y_range_tree.bit_vec.bit_vec[key-10:key+10])
                assert False
            cluster.join(cluster_dict[key])
            #print("del", key, cluster_dict[key])
            del cluster_dict[key]
        new_key = max(y_range_tree.set_to(sv_jmp.to_start(), sv_jmp.to_end(), 1, sv_jmp.from_start_same_strand()) - 1, 0)
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
        key = y_range_tree.find_zero(sv_jmp.to_start(), sv_jmp.from_start_same_strand())
        if key not in cluster_dict:
            print("CRITICAL:", key, "not in dict")
            print("searched from", y_range_tree.to_new_coord_system_c(
                sv_jmp.to_start(), sv_jmp.from_start_same_strand()))
            print(
                "bitvec:", y_range_tree.bit_vec.bit_vec[key-10:key+10])
            assert False
        cluster_dict[key].count -= 1
        if len(cluster_dict[key]) <= 0:
            # check for acceptance:
            # @note these parameters are hardcoded in two locations @todo
            if cluster_dict[key].score >= 0.3 and len(cluster_dict[key].call.supporing_jump_ids) >= 5:
                for accepted_cluster in sweep_sv_call(cluster_dict[key]):
                    print("accepting", str(accepted_cluster))
                    call_inserter.insert_call(accepted_cluster.call)
            y_range_tree.clear_downwards(key + 1)
            #print("del", key, cluster_dict[key])
            del cluster_dict[key]

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


def sv_jumps_to_dict(sv_db, run_id):
    forw_boxes_data = []
    sw_boxes_data = []
    accepted_lines_data = []
    plus_data = []
    patch_data = []
    
    sweeper = SortedSvJumpFromSql(sv_db, run_id)
    while sweeper.has_next_start():
        jump = sweeper.get_next_start()
        x = [jump.from_start_same_strand() - 0.5,
             jump.to_start() - 0.5,
             jump.from_size() + 1,
             jump.to_size() + 1,
             jump.score(),
             str(jump.id) + " " + str(jump.score())]
        if jump.does_switch_strand():
            sw_boxes_data.append(x)
        else:
            forw_boxes_data.append(x)
        if not jump.from_fuzziness_is_rightwards():
            if not jump.to_fuzziness_is_downwards():
                patch_data.append(
                    [[jump.from_pos - 2.5, jump.from_pos + .5, jump.from_pos + .5],
                     [jump.to_pos - .5, jump.to_pos + 2.5, jump.to_pos - .5]])
            else:
                patch_data.append(
                    [[jump.from_pos - 2.5, jump.from_pos + .5, jump.from_pos + .5],
                     [jump.to_pos + .5, jump.to_pos - 2.5, jump.to_pos + .5]])
        else:
            if not jump.to_fuzziness_is_downwards():
                patch_data.append(
                    [[jump.from_pos + 2.5, jump.from_pos - .5, jump.from_pos - .5],
                     [jump.to_pos - .5, jump.to_pos + 2.5, jump.to_pos - .5]])
            else:
                patch_data.append(
                    [[jump.from_pos + 2.5, jump.from_pos - .5, jump.from_pos - .5],
                     [jump.to_pos + .5, jump.to_pos - 2.5, jump.to_pos + .5]])

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
                        "type": "line",
                        "color": "green",
                        "group": "accepted_jumps",
                        "data": accepted_lines_data
                    },
                    {
                        "type": "line",
                        "color": "black",
                        "group": "diagonal",
                        "data": [[0, 0, 3*10**9, 3*10**9]]
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

    runs_from_db = SvCallerRunsFromDb(sv_db)
    cnt = 0
    while not runs_from_db.eof():
        name = runs_from_db.name()
        calls_from_db = SvCallsFromDb(sv_db, runs_from_db.id())
        accepted_boxes_data = []
        accepted_plus_data = []
        while calls_from_db.hasNext():
            jump = calls_from_db.next()
            if jump.from_size == 1 and jump.to_size == 1:
                accepted_plus_data.append([jump.from_start,
                                            jump.to_start,
                                            name + " " + str(jump.score)])
            else:
                accepted_boxes_data.append([jump.from_start - 0.5,
                                            jump.to_start - 0.5,
                                            jump.from_size + 1,
                                            jump.to_size + 1,
                                            0,
                                            name + " " + str(jump.score) + " (" +
                                            str(len(jump.supporing_jump_ids)) + ")"])
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
        runs_from_db.next()

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
