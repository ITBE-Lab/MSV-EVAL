import random
import json

##
# @brief apply modifications
# @details
# revcomp, deletions, mutations, insertions
# the order is important.
# code makes sure that the same nucleotide is never modified twice
# also indels must be at least one nuc apart from each other...
#
def disfigure(q, prob_modifier=0.1, del_prob=0.0277, mut_prob=0.0382, ins_prob=0.0651, ins_sizes=[1], 
              del_sizes=[1]):
    #deletions
    right_shifts = []
    pos = 0
    while pos < len(q):
        if random.random() < del_prob * prob_modifier:
            d_size = random.choice(del_sizes)
            right_shifts.append((pos, d_size))
            q = q[:pos] + q[pos + d_size:]
            pos += 1
        pos += 1

    # mutations
    pos = 0
    while pos < len(q):
        if random.random() < mut_prob * prob_modifier:
            q = q[:pos] + random.choice([c for c in "acgt" if c != q[pos].lower()]) + q[pos+1:]
        pos += 1

    # insertions
    pos = 0
    up_shifts = []
    while pos < len(q):
        if random.random() < ins_prob * prob_modifier:
            ins_size = random.choice(ins_sizes)
            up_shifts.append((pos,ins_size))
            ins_str = ""
            for _ in range(ins_size):
                ins_str += random.choice("acgt")
            q = q[:pos] + ins_str + q[pos:]
            pos += ins_size + 1
        else:
            pos += 1
    return q, up_shifts, right_shifts

def load_sampled_from_file(file):
    def _decode(o):
        if isinstance(o, str) or isinstance(o, unicode):
            try:
                return float(o)
            except ValueError:
                return o
        elif isinstance(o, dict):
            return {_decode(k): _decode(v) for k, v in o.items()}
        elif isinstance(o, list):
            return [_decode(v) for v in o]
        else:
            return o
    with open(file, "r") as f:
        json_file = json.loads(f.read(), object_hook=_decode)
        ins_prob, del_prob, mut_prob, ins_length_distrib, del_length_distrib, read_length_distrib = json_file
    return ins_prob, del_prob, mut_prob, ins_length_distrib, del_length_distrib, read_length_distrib