from sv_jump import *
from MA import *
import sqlite3
import math


def compute_sv_jumps(parameter_set_manager, conn, fm_index, sv_db):
    cur = conn.cursor()
    seeding_module = BinarySeeding(parameter_set_manager)

    num_destinations = parameter_set_manager.get_selected().by_name(
        "num destinations").get()
    ## fuzziness = parameter_set_manager.get_selected().by_name("fuzziness").get()

    sv_jumps = []
    jump_inserter = SvJumpInserter(sv_db, "MA-SV", "python implementation of MA-SV")

    # iterate over all queries
    cur.execute("SELECT read_table.sequence, read_table.id FROM read_table")
    for nuc_seq_blob, nuc_seq_id in cur:
        # extract seeds for query
        nuc_seq = nuc_seq_from_bytes(nuc_seq_blob)
        read_context = jump_inserter.insert_read(nuc_seq)
        segments = seeding_module.execute(fm_index, nuc_seq)
        # extract seeds & convert to python list
        seeds = [x for x in segments.extract_seeds(
            fm_index, 100, 18, len(nuc_seq), True)]
        # sort by query positions
        seeds.sort(key=lambda x: x.start)

        # compute sv jumps for all seeds
        for i, seed in enumerate(seeds):
            # create the jump from the start of the seed (temporary generator object)
            from_start_jump = PySvJump(seed, nuc_seq_id, len(nuc_seq), True)
            # as well as the one from the end of the seed (temporary generator object)
            from_end_jump = PySvJump(seed, nuc_seq_id, len(nuc_seq), False)

            # fill in at most "num_destinations" destination seeds (previous seeds on query)
            if i >= 1:
                jumps = [num_destinations, num_destinations]
                for dest_seed in seeds[i-1::-1]:
                    if max(jumps) <= 0:
                        break
                    if dest_seed.size >= 30:
                        jumps[1] -= 1
                    elif jumps[0] <= 0:
                        continue
                    jumps[0] -= 1
                    d = from_start_jump.add_destination(dest_seed)
                    if not d is None:
                        #if nuc_seq_id == 2:
                        #    print(d)
                        sv_jumps.append(d)
                        sv_jump_cpp = SvJump(seed, dest_seed, True)
                        if sv_jump_cpp.from_pos != d.x:
                            print(sv_jump_cpp.from_start, "!=", d.x)
                            assert False
                        if sv_jump_cpp.to_pos != d.y:
                            print(sv_jump_cpp.from_start, "!=", d.y)
                            assert False
                        if sv_jump_cpp.from_fuzziness_is_rightwards() != (d.fuzziness_from_dir == "right"):
                            print(sv_jump_cpp.from_fuzziness_is_rightwards(), "!=", d.fuzziness_from_dir)
                            assert False
                        if sv_jump_cpp.to_fuzziness_is_downwards() != (d.fuzziness_to_dir == "down"):
                            print(sv_jump_cpp.to_fuzziness_is_downwards(), "!=", d.fuzziness_to_dir)
                            assert False
                        read_context.insert_jump(sv_jump_cpp)
            # fill in at most "num_destinations" destination seeds (next seeds on query)
            jumps = [num_destinations, num_destinations]
            for dest_seed in seeds[i+1:]:
                if max(jumps) <= 0:
                    break
                if dest_seed.size >= 30:
                    jumps[1] -= 1
                elif jumps[0] <= 0:
                    continue
                jumps[0] -= 1
                d = from_end_jump.add_destination(dest_seed)
                if not d is None:
                    #if nuc_seq_id == 2:
                    #    print(d)
                    sv_jumps.append(d)
                    sv_jump_cpp = SvJump(seed, dest_seed, False)
                    if sv_jump_cpp.from_pos != d.x:
                        print(sv_jump_cpp.from_start, "!=", d.x)
                        assert False
                    if sv_jump_cpp.to_pos != d.y:
                        print(sv_jump_cpp.from_start, "!=", d.y)
                        assert False
                    if sv_jump_cpp.from_fuzziness_is_rightwards() != (d.fuzziness_from_dir == "right"):
                        print(sv_jump_cpp.from_fuzziness_is_rightwards(), "!=", d.fuzziness_from_dir)
                        assert False
                    if sv_jump_cpp.to_fuzziness_is_downwards() != (d.fuzziness_to_dir == "down"):
                        print(sv_jump_cpp.to_fuzziness_is_downwards(), "!=", d.fuzziness_to_dir)
                        assert False
                    read_context.insert_jump(sv_jump_cpp)
    return sv_jumps
