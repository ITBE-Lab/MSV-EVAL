from sv_jump import *
from MA import *
import sqlite3
import math


def compute_sv_jumps(parameter_set_manager, conn, fm_index, pack, sv_db):
    cur = conn.cursor()

    num_destinations = parameter_set_manager.get_selected().by_name(
        "num destinations").get()
    ## fuzziness = parameter_set_manager.get_selected().by_name("fuzziness").get()

    sv_jumps = []

    # iterate over all queries
    cur.execute(""" SELECT read_table.sequence, read_table.name 
                    FROM read_table 
                    WHERE read_table.id NOT IN ( 
                       SELECT paired_read_table.first_read FROM paired_read_table  
                       UNION SELECT paired_read_table.second_read FROM paired_read_table
                    ) """)
    nucSeqVec = ContainerVectorNucSeq()
    for nuc_seq_blob, name in cur:
        # extract seeds for query
        nuc_seq = nuc_seq_from_bytes(nuc_seq_blob)
        nuc_seq.name = name
        nucSeqVec.append(nuc_seq)

    conn.close()

    splitter_module = NucSeqSplitter(parameter_set_manager)
    lock_module = Lock(parameter_set_manager)
    seeding_module = BinarySeeding(parameter_set_manager)
    jumps_from_seeds = SvJumpsFromSeeds(parameter_set_manager)
    jumps_to_db = SvDbInserter(parameter_set_manager, sv_db, "python built compt graph")

    nucSeqVec_pledge = Pledge()
    nucSeqVec_pledge.set(nucSeqVec)
    fm_pledge = Pledge()
    fm_pledge.set(fm_index)
    pack_pledge = Pledge()
    pack_pledge.set(pack)

    
    res = VectorPledge()
    for _ in range(parameter_set_manager.get_num_threads()):
        queries_pledge = promise_me(lock_module, promise_me(splitter_module, nucSeqVec_pledge))
        unlock_module = UnLock(parameter_set_manager, queries_pledge)
        segments_pledge = promise_me(seeding_module, fm_pledge, queries_pledge)
        jumps_pledge = promise_me(jumps_from_seeds, segments_pledge, pack_pledge, fm_pledge)
        write_to_db_pledge = promise_me(jumps_to_db, jumps_pledge, queries_pledge)
        unlock_pledge = promise_me(unlock_module, write_to_db_pledge)
        res.append(unlock_pledge)
    res.simultaneous_get(parameter_set_manager.get_num_threads())

    exit(0)

