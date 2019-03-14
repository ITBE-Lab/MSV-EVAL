from sql_interface import *

def compute_socs(sv_db, pack, fm_index, parameter_set_manager=ParameterSetManager(), force=True):
    if force:
        print("forcing soc re-computation")
        sv_db.clear_soc_table()
    if sv_db.has_socs():
        print("found SoC's -> will not recompute them")
        return

    fm_pledge = Pledge()
    fm_pledge.set(fm_index)
    pack_pledge = Pledge()
    pack_pledge.set(pack)
    read_pledge = sv_db.get_reads_pledge(parameter_set_manager)

    module_lock = Lock(parameter_set_manager)
    seeding_module = BinarySeeding(parameter_set_manager)
    soc_module = StripOfConsideration(parameter_set_manager)
    soc_inserter_module = sv_db.get_soc_inserter_module(parameter_set_manager)

    res_vec = VectorPledge()
    for _ in range(parameter_set_manager.get_num_threads()):
        query = promise_me(module_lock, read_pledge)
        seeds = promise_me(seeding_module, fm_pledge, query)
        soc_queue = promise_me(soc_module, seeds, query, pack_pledge, fm_pledge)
        empty = promise_me(soc_inserter_module, query, soc_queue)

        res_vec.append(promise_me(UnLock(parameter_set_manager, query), empty))

    res_vec.simultaneous_get(parameter_set_manager.get_num_threads())



if __name__ == "__main__":
    fm_index = FMIndex()
    fm_index.load("/MAdata/genome/human/GRCh38.p12/ma/genome")
    pack = Pack()
    pack.load("/MAdata/genome/human/GRCh38.p12/ma/genome")
    compute_socs(SV_DB("/MAdata/databases/sv_simulated"), pack, fm_index, force=True)