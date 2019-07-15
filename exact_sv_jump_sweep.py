from MA import *
from svCallPy import *
from bokeh.plotting import figure, show, reset_output, ColumnDataSource
PRINTS = False

def sweep_sv_jumps(sv_jmps, estimated_coverage, re_estimate_cluster_size=True):
    def to_end(sv_jmp):
        return sv_jmp.to_end() - 1 if sv_jmp.switch_strand_known() else sv_jmp.to_start() + sv_jmp.from_size()

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
        if len(cluster.call.supporing_jump_ids) >= max(2, estimated_coverage/8):
            if re_estimate_cluster_size:
                right = cluster.right()
                up = cluster.up()
                cluster.call.from_start = max(0, cluster.left())
                cluster.call.to_start = max(0, cluster.down())
                cluster.call.from_size = max(right - cluster.call.from_start, 1)
                cluster.call.to_size = max(1, up - cluster.call.to_start)
            cluster.call.num_supp_nt = 0
            for x in range(len(cluster.call.supporing_jump_ids)):
                cluster.call.num_supp_nt += cluster.call.get_jump(x).num_supp_nt()
            # if the cluster still fulfills the required criteria
            # @note another hardcoded quality parameter @todo
            #if cluster.call.score >= 500:
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

## complete linkage clustering for jump distances
# we call sweep_sv_jumps for all insert_ratio clusters with a max dist of max_insert_ratio_diff
# this clustering is necessary because there might be an edge in the graph that has several different inserted
# sequences. We need to consider these sequences individually -> cluster by sequence length
# if the sequences are different by nucleotides, we need to figure that out later in the multialignment step...
#
def sweep_sv_call(sv_call, estimated_coverage):
    max_insert_ratio_diff = 150

    sv_jmps = [sv_call.call.get_jump(x) for x in range(len(sv_call.call.supporing_jump_ids))]
    sv_jmps.sort(key=lambda x: (x.insert_ratio(), x.query_distance()))

    # i & j are the start and end positions of the clusters, respectiveley
    # this works by setting i to the start of the cluster and then gradually increasing j while it belongs to the 
    # complete linkage cluster (we work on sorted linear data)
    i = 0
    j = 0
    ret = []
    while i < len(sv_jmps):
        # increase j if the insert_ratio between the 'i' and 'j' object is closer than 'max_insert_ratio_diff'
        # if we reach a tail edge (those edges are sorted to the end since their insert 
        # size is 'inf') we check if the current insert ratio is larger than the tail of the read that created the edge
        # if so we join the tail edge into the cluster.
        if j < len(sv_jmps) and sv_jmps[i].insert_ratio() >= \
                    (sv_jmps[j].insert_ratio() - max_insert_ratio_diff if sv_jmps[j].switch_strand_known() else sv_jmps[j].query_distance()):
            j += 1
        else:
            # perform the overlapping on the cluster
            ret.extend(sweep_sv_jumps(sv_jmps[i:j], estimated_coverage))
            # set i to the start of the next cluster
            i = j
    return ret
