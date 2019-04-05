from sv_jump import *
from MA import *
import sqlite3
import math


def compute_sv_jumps(parameter_set_manager, conn, fm_index):
    cur = conn.cursor()
    seeding_module = BinarySeeding(parameter_set_manager)

    num_destinations = parameter_set_manager.get_selected().by_name(
        "num destinations").get()
    ## fuzziness = parameter_set_manager.get_selected().by_name("fuzziness").get()

    sv_jumps = []

    # iterate over all queries
    cur.execute("SELECT read_table.sequence, read_table.id FROM read_table")
    for nuc_seq_blob, nuc_seq_id in cur:
        # extract seeds for query
        nuc_seq = nuc_seq_from_bytes(nuc_seq_blob)
        segments = seeding_module.execute(fm_index, nuc_seq)
        # extract seeds & convert to python list
        seeds = [x for x in segments.extract_seeds(
            fm_index, 100, 19, len(nuc_seq), False)]
        # sort by query positions
        seeds.sort(key=lambda x: x.start)
        # if nuc_seq_id == 13:
        #    for seed in seeds:
        #        print(seed.start, seed.start + seed.size, seed.start_ref, seed.start_ref + seed.size, seed.on_forward_strand)

        # compute sv jumps for all seeds
        for i, seed in enumerate(seeds):
            # create the jump from the start of the seed (temporary generator object)
            from_start_jump = SvJump(seed, nuc_seq_id, len(nuc_seq), True)
            # as well as the one from the end of the seed (temporary generator object)
            from_end_jump = SvJump(seed, nuc_seq_id, len(nuc_seq), False)

            # fill in at most "num_destinations" destination seeds (previous seeds on query)
            if i >= 1:
                for dest_seed in seeds[max(i-num_destinations-1, 0):i]:
                    d = from_start_jump.add_destination(dest_seed)
                    if not d is None:
                        sv_jumps.append(d)
            # fill in at most "num_destinations" destination seeds (next seeds on query)
            for dest_seed in seeds[i+1:i+num_destinations+1]:
                d = from_end_jump.add_destination(dest_seed)
                if not d is None:
                    sv_jumps.append(d)
    return sv_jumps
