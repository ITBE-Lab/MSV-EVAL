from MA import *
from MSV import *


def comp_mems(query_genome, reference_genome):
    # load reference from file    
    ref_pack = Pack()
    ref_pack.load(reference_genome + "/ma/genome")

    # load query from file    
    q_pack = Pack()
    q_pack.load(query_genome + "/ma/genome")

    # cut out region
    # see bottom of libs/ma/src/container/pack.cpp for functions of Pack()
    query = q_pack.extract_from_to(10000, 10100)

    # set parameters
    param = ParameterSetManager()
    param.by_name("Minimizers - k").set(15)
    param.by_name("Minimizers - w").set(10)

    # compute mm index
    mm_index = MinimizerIndex(param, libMA.util.StringVector([str(ref_pack.extract_forward_strand())]),
                              libMA.util.StringVector(["fullGenome"]))

    # create module objects
    # get minimizers
    seeding_module = MinimizerSeeding(param)
    # turn minimizers into MEMs
    seed_lumper = SeedLumping(param)

    # get mems
    if False:
        minimizers = seeding_module.execute(mm_index, query, ref_pack)
        mems = seed_lumper.execute(minimizers, query, ref_pack)
    else:
        mm_index_pledge = Pledge()
        mm_index_pledge.set(mm_index)

        ref_pack_pledge = Pledge()
        ref_pack_pledge.set(ref_pack)

        query_pledge = Pledge()
        query_pledge.set(query)

        minimizer_pledge = promise_me(seeding_module, mm_index_pledge, query_pledge, ref_pack_pledge)
        mems_pledge = promise_me(seed_lumper, minimizer_pledge, query_pledge, ref_pack_pledge)

        mems = mems_pledge.get()

        
        #query_pledge.set(other_query)
        #other_mems = mems_pledge.get()


    # iterate over mems
    for mem in mems:
        print("q,r,l:", mem.start, mem.start_ref, mem.size)


comp_mems("/MAdata/genome/yeasts/UFRJ50816", "/MAdata/genome/yeasts/YPS138")
