from sv_jump import *
from MA import *
import sqlite3

def compute_sv_jumps(parameter_set_manager, conn, fm_index):
    cur = conn.cursor()
    seeding_module = BinarySeeding(parameter_set_manager)

    num_destinations = parameter_set_manager.get_selected().by_name("num destinations").get()
    ## fuzziness = parameter_set_manager.get_selected().by_name("fuzziness").get()

    sv_jumps = []

    # iterate over all queries
    cur.execute("SELECT read_table.sequence, read_table.id FROM read_table")
    for nuc_seq_blob, nuc_seq_id in cur:
        # extract seeds for query
        nuc_seq = nuc_seq_from_bytes(nuc_seq_blob)
        segments = seeding_module.execute(fm_index, nuc_seq)
        # extract seeds & convert to python list
        seeds = [x for x in segments.extract_seeds(fm_index, 100, 0, len(nuc_seq), False)]
        # sort by query positions
        seeds.sort(key=lambda x: x.start)
        #if nuc_seq_id == 13:
        #    for seed in seeds:
        #        print(seed.start, seed.start + seed.size, seed.start_ref, seed.start_ref + seed.size, seed.on_forward_strand)

        # compute sv jumps for all seeds
        for i, seed in enumerate(seeds):
            # add the jump itself
            sv_jumps.append(SvJump(seed, nuc_seq_id, len(nuc_seq)))
            # fill in at most "num_destinations" destination seeds (next seeds on query)
            rem_destination = num_destinations
            for dest_seed in seeds[i:] if seed.on_forward_strand else seeds[i::-1]:
                if sv_jumps[-1].add_destination(dest_seed): # function returns wether seed was actually added
                    rem_destination -= 1
                    # if we have added enough destination seeds, break the loop
                    if rem_destination == 0:
                        break
    return sv_jumps
