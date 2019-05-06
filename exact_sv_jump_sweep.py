from MA import *
from svCallPy import *

def sweep_sv_jumps(sv_jmps, re_estimate_cluster_size=True):
    if len(sv_jmps) == 0: # sanity check
        return
    # squash sv_jump indices
    pos_list = []
    for sv_jmp in sv_jmps:
        pos_list.append(sv_jmp.to_start())
        pos_list.append(sv_jmp.to_end())
    pos_list.sort()
    pos_dict = {}
    idx = 0
    for pos in pos_list:
        if not pos in pos_dict:
            pos_dict[pos] = idx
            idx += 1

    # create start list
    sv_jmps_start = []
    for sv_jmp in sv_jmps:
        sv_jmps_start.append(sv_jmp)
    sv_jmps_start.sort(key=lambda x: x.from_start())

    # create end list
    sv_jmps_end = []
    for sv_jmp in sv_jmps:
        sv_jmps_end.append(sv_jmp)
    sv_jmps_end.sort(key=lambda x: x.from_end())

    # setup sweep vector
    sweep_vec = []
    for _ in range(idx+1):
        sweep_vec.append([0,None])

    # define linesweep helpers
    def helper_start(sv_jmp):
        cluster = SvCallPy(sv_jmp)
        # join with all overlapping clusters
        i = pos_dict[sv_jmp.to_start()]
        while i <= pos_dict[sv_jmp.to_end()]:
            if sweep_vec[i][0] > 0:
                cluster.join(sweep_vec[i][1])
                # we can jump to the end of the discovered cluster immediately
                i = max(i, pos_dict[sweep_vec[i][1].call.to_size + sweep_vec[i][1].call.to_start - 1])
            i += 1
        # insert the newly computed cluster above and below the current sv_jump (we can't know the actual cluster size)
        i = pos_dict[cluster.call.to_start]
        while i < pos_dict[sv_jmp.to_start()]:
            if sweep_vec[i][1] == sweep_vec[pos_dict[sv_jmp.to_start()]][1]:
                sweep_vec[i][1] = cluster
            i += 1
        i = pos_dict[sv_jmp.to_end()] + 1
        while i <= pos_dict[cluster.call.to_size + cluster.call.to_start - 1]:
            if sweep_vec[i][1] == sweep_vec[pos_dict[sv_jmp.to_end()]][1]:
                sweep_vec[i][1] = cluster
            i += 1

        # increment the sv_jump counter and insert the cluster overlapping with the current
        i = pos_dict[sv_jmp.to_start()]
        while i <= pos_dict[sv_jmp.to_end()]:
            sweep_vec[i][0] += 1
            sweep_vec[i][1] = cluster
            i += 1

    ret = []
    def check_cluster(cluster):
        # remove jumps with equal read id's
        new_jmps = {}
        for x in range(len(cluster.call.supporing_jump_ids)):
            jmp = cluster.call.get_jump(x)
            if not jmp.read_id in new_jmps:
                new_jmps[jmp.read_id] = jmp
            elif jmp.query_distance() < new_jmps[jmp.read_id].query_distance():
                new_jmps[jmp.read_id] = jmp

        cluster.call.clear_jumps()
        for jmp in new_jmps.values():
            cluster.call.add_jump(jmp)

        # if the cluster still fulfills the required criteria
        # @note these parameters are hardcoded in two locations @todo
        if cluster.score >= 0.3 and len(cluster.call.supporing_jump_ids) >= 5:
            if re_estimate_cluster_size:
                right = cluster.right()
                up = cluster.up()
                cluster.call.from_start = cluster.left()
                cluster.call.to_start = cluster.down()
                cluster.call.from_size = right - cluster.call.from_start
                cluster.call.to_size = up - cluster.call.to_start
            ret.append(cluster)

    def helper_end(sv_jmp):
        # get one pointer to the current cluster...
        i = pos_dict[sv_jmp.to_start()]
        assert sweep_vec[i][0] > 0
        # decrease the amount ef elements in the cluster
        sweep_vec[i][1].count -= 1
        # if the count hits zero check if the cluster is worth keeping
        if len(sweep_vec[i][1]) == 0:
            try:
                check_cluster(sweep_vec[i][1])
            except:
                print(pos_dict.items(), i)
                print("x", sv_jmp.from_start(), sv_jmp.from_end(), sv_jmp.to_start(), sv_jmp.to_end())
                for sv_jmp in sv_jmps:
                    print(sv_jmp.from_start(), sv_jmp.from_end(), sv_jmp.to_start(), sv_jmp.to_end())
                for idx, sv_jmp in sweep_vec:
                    print(idx, sv_jmp)
                assert False
        # decrement the sv_jump counter
        i = pos_dict[sv_jmp.to_start()]
        while i <= pos_dict[sv_jmp.to_end()]:
            sweep_vec[i][0] -= 1
            i += 1

    # do the actual sweep:
    i = 0
    j = 0
    while i < len(sv_jmps_start) and j < len(sv_jmps_end):
        if sv_jmps_start[i].from_start() <= sv_jmps_end[j].from_end():
            helper_start(sv_jmps_start[i])
            i += 1
        else:
            helper_end(sv_jmps_end[j])
            j += 1
    assert i == len(sv_jmps_start)
    while j < len(sv_jmps_end):
        helper_end(sv_jmps_end[j])
        j += 1

    for i, x in sweep_vec:
        assert i == 0

    return ret

def sweep_sv_call(sv_call):
    # single linkage clustering for jump distances
    # we call sweep_sv_jumps for all insert_ratio clusters with a max dist of max_insert_ratio_diff
    max_insert_ratio_diff = 500

    sv_jmps = [sv_call.call.get_jump(x) for x in range(len(sv_call.call.supporing_jump_ids))]
    sv_jmps.sort(key=lambda x: x.insert_ratio())

    i = 0
    j = 0
    while j < len(sv_jmps):
        if sv_jmps[i].insert_ratio() + max_insert_ratio_diff >= sv_jmps[j].insert_ratio():
            j += 1
        else:
            # @note currently we always accept the sv_jump with the smaller insert size
            # otherwise we would create artifacts (due to the multiple seeds beeing used)...
            # @todo is there a smarter way to filter the artifacts 
            #       or is it possible to have two insertions with different lengths at the same spot?
            ret = sweep_sv_jumps(sv_jmps[i:j])
            if len(ret) > 0:
                return ret
            i = j
    return sweep_sv_jumps(sv_jmps[i:j])
