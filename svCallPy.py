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
        if jump.from_fuzziness_is_rightwards():
            # if fuzziness is rightswards, this jump indicates the left border...
            self.l_left.append( jump.from_pos )
        else:
            self.l_right.append( jump.from_pos )
        if jump.to_fuzziness_is_downwards():
            self.l_up.append( jump.to_pos )
        else:
            self.l_down.append( jump.to_pos )
        self.t = 5 # confidence in clustersize optimization

    def right(self):
        if len(self.l_right) < 3:
            return self.call.from_size + self.call.from_start
        self.l_right.sort()
        return self.l_right[int(len(self.l_right)/self.t)] + 5

    def left(self):
        if len(self.l_left) < 3:
            return self.call.from_start
        self.l_left.sort(reverse=True)
        return self.l_left[int(len(self.l_left)/self.t)] - 5

    def up(self):
        if len(self.l_left) < 3:
            return self.call.to_size + self.call.to_start
        self.l_up.sort()
        return self.l_up[int(len(self.l_up)/self.t)] + 5

    def down(self):
        if len(self.l_down) < 3:
            return self.call.to_start
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