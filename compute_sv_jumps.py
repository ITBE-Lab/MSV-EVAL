from sv_jump import *
from MA import *
import sqlite3
import math

def compute_sv_jumps(parameter_set_manager, conn, fm_index, pack, sv_db):
    # iterate over all queries
    cur = conn.cursor()
    cur.execute("""SELECT read_table.sequence, read_table.id FROM read_table""")
    nucSeqVec = ContainerVectorNucSeq()
    for nuc_seq_blob, nuc_seq_id in cur:
        # extract seeds for query
        nuc_seq = nuc_seq_from_bytes(nuc_seq_blob)
        nuc_seq.name = "read_nr_" + str(nuc_seq_id)
        nucSeqVec.append(nuc_seq)

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
        segments_pledge = promise_me(seeding_module, fm_pledge, queries_pledge)
        jumps_pledge = promise_me(jumps_from_seeds, segments_pledge, pack_pledge, fm_pledge)
        write_to_db_pledge = promise_me(jumps_to_db, jumps_pledge, queries_pledge)
        unlock_pledge = promise_me(UnLock(parameter_set_manager, queries_pledge), write_to_db_pledge)
        res.append(unlock_pledge)
    res.simultaneous_get(1 )#parameter_set_manager.get_num_threads()) @todo check for parallel bug...

