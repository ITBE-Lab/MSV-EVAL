from MA import *
from svCallPy import *
from bokeh.plotting import figure, show, reset_output, ColumnDataSource
PRINTS = False

def sweep_sv_jumps(sv_jmps, re_estimate_cluster_size=True):
    def to_end(sv_jmp):
        return sv_jmp.to_end() if sv_jmp.switch_strand_known() else sv_jmp.to_start() + sv_jmp.from_size()

    if len(sv_jmps) == 0: # sanity check
        return
    # squash sv_jump indices
    pos_list = []
    if PRINTS:
        print("==================")
    for sv_jmp in sv_jmps:
        if PRINTS:
            print(sv_jmp.from_start(), sv_jmp.to_start(), sv_jmp.from_end(), sv_jmp.to_end())
        pos_list.append(sv_jmp.to_start())
        pos_list.append(to_end(sv_jmp))
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
    for _ in range(idx):
        sweep_vec.append([0,None])

    # define linesweep helpers
    def helper_start(sv_jmp):
        cluster = SvCallPy(sv_jmp)
        if not sv_jmp.switch_strand_known():
            cluster.call.to_size = sv_jmp.from_size() + 1
        # join with all overlapping clusters
        i = pos_dict[sv_jmp.to_start()]
        joined_clusters = set()
        while i <= pos_dict[to_end(sv_jmp)]:
            if sweep_vec[i][0] > 0 and not sweep_vec[i][1] in joined_clusters:
                cluster.join(sweep_vec[i][1])
                joined_clusters.add(sweep_vec[i][1])
                # we can jump to the end of the discovered cluster immediately
                # NO WE CANNOT -> the cluster may have a u shape... like this: 
                #  xxxxx
                #  x
                #  xxxxx
                #i = max(i, pos_dict[sweep_vec[i][1].call.to_size + sweep_vec[i][1].call.to_start - 1])
            i += 1
        # insert the newly computed cluster above and below the current sv_jump (we can't know the actual cluster size)
        i = pos_dict[cluster.call.to_start]
        while i < pos_dict[sv_jmp.to_start()]:
            if sweep_vec[i][1] in joined_clusters:
                sweep_vec[i][1] = cluster
            i += 1
        i = pos_dict[to_end(sv_jmp)] + 1
        while i <= pos_dict[cluster.call.to_size + cluster.call.to_start - 1]:
            if sweep_vec[i][1] in joined_clusters:
                sweep_vec[i][1] = cluster
            i += 1

        # increment the sv_jump counter and insert the cluster overlapping with the current
        i = pos_dict[sv_jmp.to_start()]
        while i <= pos_dict[to_end(sv_jmp)]:
            sweep_vec[i][0] += 1
            sweep_vec[i][1] = cluster
            i += 1

        if PRINTS:
            print( *(x for x,y in sweep_vec), ":: start ::", i)
            print( *(y.count if not y is None else 0 for x,y in sweep_vec), ":: cluster ::")

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
        if cluster.score >= 5 and len(cluster.call.supporing_jump_ids) >= 5:
            if re_estimate_cluster_size:
                right = cluster.right()
                up = cluster.up()
                cluster.call.from_start = cluster.left()
                cluster.call.to_start = cluster.down()
                cluster.call.from_size = max(right - cluster.call.from_start, 1)
                cluster.call.to_size = max(1, up - cluster.call.to_start)
            ret.append(cluster)

    def helper_end(sv_jmp):
        # get one pointer to the current cluster...
        i = pos_dict[sv_jmp.to_start()]
        assert sweep_vec[i][0] > 0
        # decrease the amount ef elements in the cluster
        sweep_vec[i][1].count -= 1
        # if the count hits zero check if the cluster is worth keeping
        try:
            if len(sweep_vec[i][1]) == 0:
                check_cluster(sweep_vec[i][1])
        except:
            plot = figure()
            for sv_jmp_ in sv_jmps:
                plot.quad(sv_jmp_.from_start(), sv_jmp_.from_end(), sv_jmp_.to_start(), sv_jmp_.to_end(), fill_alpha=0.5)
            plot.quad(sv_jmp.from_start(), sv_jmp.from_end(), sv_jmp.to_start(), sv_jmp.to_end(),
                      fill_alpha=0.5, color="red")
            show(plot)
            #print(pos_dict.items(), i)
            #print("x", sv_jmp.from_start(), sv_jmp.from_end(), sv_jmp.to_start(), to_end(sv_jmp))
            #for sv_jmp in sv_jmps:
            #    print(sv_jmp.from_start(), sv_jmp.from_end(), sv_jmp.to_start(), to_end(sv_jmp))
            #for idx, sv_jmp in sweep_vec:
            #    print(idx, sv_jmp)
            assert False
        # decrement the sv_jump counter
        i = pos_dict[sv_jmp.to_start()]
        while i <= pos_dict[to_end(sv_jmp)]:
            sweep_vec[i][0] -= 1
            i += 1

        if PRINTS:
            print( *(x for x,y in sweep_vec), ":: end ::", i)
            print( *(y.count if not y is None else 0 for x,y in sweep_vec), ":: cluster ::")

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

#
# @todo problem: illumina reads have a insert ratio '> remaining query distance'
# this should mean that they get added to each cluster that matches that
# -> implement that
def sweep_sv_call(sv_call):
    # single linkage clustering for jump distances
    # we call sweep_sv_jumps for all insert_ratio clusters with a max dist of max_insert_ratio_diff
    max_insert_ratio_diff = 500
    max_q_len = 200

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
