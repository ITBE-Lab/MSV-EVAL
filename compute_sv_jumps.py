from sv_jump import *
from MA import *
import sqlite3
import math
from analyze_runtimes import AnalyzeRuntimes

def compute_sv_jumps(parameter_set_manager, fm_index, pack, sv_db, seq_id=0):
    lock_module = Lock(parameter_set_manager)
    seeding_module = BinarySeeding(parameter_set_manager)
    jumps_to_db = SvDbInserter(parameter_set_manager, sv_db, "python built compt graph")
    jumps_from_seeds = SvJumpsFromSeeds(parameter_set_manager, jumps_to_db.cpp_module.jump_inserter.sv_jump_run_id,
                                        sv_db)

    fm_pledge = Pledge()
    fm_pledge.set(fm_index)
    pack_pledge = Pledge()
    pack_pledge.set(pack)

    analyze = AnalyzeRuntimes()

    res = VectorPledge()
    # graph for single reads
    for idx in range(parameter_set_manager.get_num_threads()):
        nuc_seq_getter = AllNucSeqFromSql(parameter_set_manager, sv_db, seq_id, idx,
                                          parameter_set_manager.get_num_threads())
        queries_pledge = promise_me(nuc_seq_getter)
        analyze.register("AllNucSeqFromSql", queries_pledge)
        query_pledge = promise_me(lock_module, queries_pledge)
        segments_pledge = promise_me(seeding_module, fm_pledge, query_pledge)
        analyze.register("BinarySeeding", segments_pledge)
        jumps_pledge = promise_me(jumps_from_seeds, segments_pledge, pack_pledge, fm_pledge, query_pledge)
        analyze.register("SvJumpsFromSeeds", jumps_pledge)
        write_to_db_pledge = promise_me(jumps_to_db, jumps_pledge, query_pledge)
        analyze.register("SvDbInserter", write_to_db_pledge)
        unlock_pledge = promise_me(UnLock(parameter_set_manager, query_pledge), write_to_db_pledge)
        res.append(unlock_pledge)

    # drain all sources
    res.simultaneous_get( parameter_set_manager.get_num_threads() )

    analyze.analyze()
    
    sv_db.create_jump_indices( jumps_to_db.cpp_module.jump_inserter.sv_jump_run_id )

    # return the run_id
    return jumps_to_db.cpp_module.jump_inserter.sv_jump_run_id

