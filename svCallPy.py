from MA import *

class SvCallPy:
    def __init__(self, jump):
        self.call = SvCall(jump)
        self.count = 1
        self.score = jump.score()
        self.l_right = []
        self.l_up = []
        self.l_down = []
        self.l_left = []
        self.del_to_ins_ratio = 1
        f = jump.from_pos
        t = jump.to_pos
        if not jump.from_known():
            f = t
        if not jump.to_known():
            t = f
        if jump.from_fuzziness_is_rightwards():
            # if fuzziness is rightswards, this jump indicates the left border...
            self.l_left.append( f )
            x = f + int(min( jump.ref_distance() / self.del_to_ins_ratio,
                                     jump.query_distance() * self.del_to_ins_ratio ))
            self.l_right.append( x )
        else:
            self.l_right.append( f )
            x = f - int(min( jump.ref_distance() / self.del_to_ins_ratio,
                                     jump.query_distance() * self.del_to_ins_ratio ))
            self.l_left.append( x )
        if jump.to_fuzziness_is_downwards():
            self.l_up.append( t )
            x = t - int(min( jump.ref_distance() / self.del_to_ins_ratio,
                                     jump.query_distance() * self.del_to_ins_ratio ))
            self.l_down.append( x )
        else:
            self.l_down.append( t )
            x = t + int(min( jump.ref_distance() / self.del_to_ins_ratio,
                                     jump.query_distance() * self.del_to_ins_ratio ))
            self.l_up.append( x )
        self.t = 5 # confidence in clustersize optimization

    def right(self):
        self.l_right.sort()
        return self.l_right[int(len(self.l_right)/self.t)] + 5

    def left(self):
        self.l_left.sort(reverse=True)
        return self.l_left[int(len(self.l_left)/self.t)] - 5

    def up(self):
        self.l_up.sort()
        return self.l_up[int(len(self.l_up)/self.t)] + 5

    def down(self):
        self.l_down.sort(reverse=True)
        return self.l_down[int(len(self.l_down)/self.t)] - 5

    def join(self, other):
        self.call.join(other.call)
        self.count += other.count
        self.score += other.score
        self.l_down.extend(other.l_down)
        self.l_up.extend(other.l_up)
        self.l_right.extend(other.l_right)
        self.l_left.extend(other.l_left)

    def __len__(self):
        return self.count

    def __str__(self):
        q_dist = 0
        for x in range(len(self.call.supporing_jump_ids)):
            q_dist += self.call.get_jump(x).query_distance()
        q_dist /= len(self.call.supporing_jump_ids)
        return "from " + str(self.call.from_start) + "~" + str(self.call.from_size + self.call.from_start) + " to " + str(self.call.to_start) + "~" + str(self.call.to_size + self.call.to_start) + (" switch_strand" if self.call.switch_strand else "") + " score " + str(self.score) + " cnt " + str(len(self.call.supporing_jump_ids))+ " avg. insert " + str(q_dist)