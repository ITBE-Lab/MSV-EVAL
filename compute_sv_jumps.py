from sv_jump import *
from MA import *
import sqlite3
import math

def compute_sv_jumps(parameter_set_manager, fm_index, pack, sv_db):
    parameter_set_manager.by_name("Mean Distance of Paired Reads").set(500) # @todo sample this...
    parameter_set_manager.by_name("Do Mate Jumps").set(True)
    sv_db.drop_caller_indices()
    nuc_seq_getter = None
    if parameter_set_manager.by_name("Do Mate Jumps").get():
        nuc_seq_getter = NucSeqFromSql(parameter_set_manager, sv_db)
    else:
        nuc_seq_getter = AllNucSeqFromSql(parameter_set_manager, sv_db)
    lock_module = Lock(parameter_set_manager)
    seeding_module = BinarySeeding(parameter_set_manager)
    jumps_from_seeds = SvJumpsFromSeeds(parameter_set_manager)
    jumps_to_db = SvDbInserter(parameter_set_manager, sv_db, "python built compt graph")

    fm_pledge = Pledge()
    fm_pledge.set(fm_index)
    pack_pledge = Pledge()
    pack_pledge.set(pack)

    res = VectorPledge()
    queries_pledge = promise_me(nuc_seq_getter) # @note this cannot be in the loop (synchronization!)
    # graph for single reads
    for _ in range(parameter_set_manager.get_num_threads()):
        query_pledge = promise_me(lock_module, queries_pledge)
        segments_pledge = promise_me(seeding_module, fm_pledge, query_pledge)
        jumps_pledge = promise_me(jumps_from_seeds, segments_pledge, pack_pledge, fm_pledge, query_pledge)
        write_to_db_pledge = promise_me(jumps_to_db, jumps_pledge, query_pledge)
        unlock_pledge = promise_me(UnLock(parameter_set_manager, query_pledge), write_to_db_pledge)
        res.append(unlock_pledge)
    
    # graph for paired reads
    if parameter_set_manager.by_name("Do Mate Jumps").get():
        paired_nuc_seq_getter = PairedNucSeqFromSql(parameter_set_manager, sv_db)
        paired_jumps_from_seeds = SvJumpsFromSeedsPaired(parameter_set_manager)
        module_get_first = GetFirstQuery(parameter_set_manager)
        module_get_second = GetSecondQuery(parameter_set_manager)
        paired_queries_pledge = promise_me(paired_nuc_seq_getter) # @note this cannot be in the loop (synchronization!)
        for _ in range(parameter_set_manager.get_num_threads()):
            paired_query_pledge = promise_me(lock_module, paired_queries_pledge)
            query_a_pledge = promise_me(module_get_first, paired_query_pledge)
            query_b_pledge = promise_me(module_get_second, paired_query_pledge)
            segments_a_pledge = promise_me(seeding_module, fm_pledge, query_a_pledge)
            segments_b_pledge = promise_me(seeding_module, fm_pledge, query_b_pledge)
            jumps_pledge = promise_me(paired_jumps_from_seeds, segments_a_pledge, segments_b_pledge, pack_pledge, 
                                    fm_pledge, query_a_pledge, query_b_pledge)
            write_to_db_pledge = promise_me(jumps_to_db, jumps_pledge, query_a_pledge) # @todo atm paired read jumps just get assigned to the first one...
            unlock_pledge = promise_me(UnLock(parameter_set_manager, paired_query_pledge), write_to_db_pledge)
            res.append(unlock_pledge)

    # drain all sources
    res.simultaneous_get( parameter_set_manager.get_num_threads() )
    sv_db.create_caller_indices()

    # return the run_id
    return jumps_to_db.cpp_module.jump_inserter.sv_caller_run_id

