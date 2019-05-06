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
        if not jump.from_fuzziness_is_rightwards():
            self.l_left.append( (jump.from_pos, jump.query_distance()) )
        else:
            self.l_right.append( (jump.from_pos, jump.query_distance()) )
        if not jump.to_fuzziness_is_downwards():
            self.l_up.append( (jump.to_pos, jump.query_distance()) )
        else:
            self.l_down.append( (jump.to_pos, jump.query_distance()) )

    def max_dist_list(self, l):
        max_d = max(y for x,y in l)
        ret = []
        for x, y in l:
            if max_d / 10 < y:
                ret.append(x)
        return ret

    def right(self):
        l_right = self.max_dist_list(self.l_right)
        l_right.sort(reverse=True)
        return l_right[int(len(l_right)/50)]

    def left(self):
        l_left = self.max_dist_list(self.l_left)
        l_left.sort()
        return l_left[int(len(l_left)/50)]

    def up(self):
        l_up = self.max_dist_list(self.l_up)
        l_up.sort(reverse=True)
        return l_up[int(len(l_up)/50)]

    def down(self):
        l_down = self.max_dist_list(self.l_down)
        l_down.sort()
        return l_down[int(len(l_down)/50)]

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