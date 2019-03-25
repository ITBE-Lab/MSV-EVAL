from sv_jump import *

def line_sweep_list(sorted_list, do_functor, key, fuzziness):
    sorted_list.sort(key=lambda x: key(x))
    start_it = 0
    end_it = 0

    while start_it < len(sorted_list):
        while end_it < len(sorted_list) and key(sorted_list[end_it]) - fuzziness <= key(sorted_list[start_it]):
            end_it += 1
        do_functor(sorted_list[start_it:end_it])
        start_it += 1

def sweep_sv_jumps(parameter_set_manager, sv_jumps):
    accepted_sv_jumps = []
    fuzziness = parameter_set_manager.get_selected().by_name("fuzziness").get()
    min_coverage = parameter_set_manager.get_selected().by_name("min coverage").get()

    def check_sv_jump(curr_jumps):
        read_id_set = set()
        for curr_jump in curr_jumps:
            read_id_set.add(curr_jump.read_id)
        if len(read_id_set) < min_coverage:
            return
        # extract all destinations
        dest_list = []
        for curr_jump in curr_jumps:
            dest_list.extend(curr_jump.destinations)
        
        # and sweep again
        def accept_functor(curr_destinations):
            dest_read_id_set = set()
            for curr_dest in curr_destinations:
                dest_read_id_set.add(curr_dest.parent.read_id)
            if len(dest_read_id_set) < min_coverage:
                return
            print("accepting sv jump:")
            print("from ref pos:", [x.ref_pos for x in curr_jumps])
            print("to ref pos:", [x.ref_pos for x in curr_destinations])
            print("q_distances:", [x.q_distance for x in curr_destinations])
            print("switch strand:", [x.switch_strands for x in curr_destinations])

        line_sweep_list(dest_list, accept_functor, lambda x: x.ref_pos, fuzziness)

    # line sweep
    line_sweep_list(sv_jumps, check_sv_jump, lambda x: x.ref_pos, fuzziness)

    return accepted_sv_jumps