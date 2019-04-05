from MA import *
import json
import sqlite3
NUM_THREADS = 32


def for_seed_sections(segments, fm_index, nuc_seq_len):
    seeds = [x for x in segments.extract_seeds(
        fm_index, 100, 0, nuc_seq_len, True)]
    seeds.sort(key=lambda x: x.start_ref)
    start_idx = 0
    end_idx = 1
    # a sections is a set of consecutive seeds where the maximal gap between two consecutive seeds is max_section_dist
    max_section_dist = 1000
    # yield all sections of seeds
    while start_idx < len(seeds):
        while end_idx < len(seeds) and seeds[end_idx-1].start_ref + max_section_dist >= seeds[end_idx].start_ref:
            end_idx += 1
        yield seeds[start_idx:end_idx]
        start_idx = end_idx
        end_idx += 1


def create_json_from_db(db_name, pack_path):
    parameter_manager = ParameterSetManager()
    parameter_manager.set_selected("PacBio")
    fm_index = FMIndex()
    fm_index.load("/MAdata/genome/human/GRCh38.p12/ma/genome")
    fm_pledge = Pledge()
    fm_pledge.set(fm_index)

    # set up all modules
    database = SV_DB(db_name, "open")
    sql = NucSeqFromSql(parameter_manager, database)
    lock = Lock(parameter_manager)
    paired_sql = PairedNucSeqFromSql(parameter_manager, database)
    seeding = BinarySeeding(parameter_manager)
    collector = libMA.NucSeqSegmentCollector(parameter_manager)
    module_get_first = GetFirstQuery(parameter_manager)
    module_get_second = GetSecondQuery(parameter_manager)
    join_module = libMA.ContainerJoin(parameter_manager)

    # set up the graph
    queries = promise_me(sql)
    paired_queries = promise_me(paired_sql)
    res_vec = VectorPledge()
    for _ in range(NUM_THREADS):
        # for single
        query = promise_me(lock, queries)
        seeds = promise_me(seeding, fm_pledge, query)
        empty = promise_me(collector, query, seeds)
        res_vec.append(promise_me(UnLock(parameter_manager, query), empty))

        # for paired
        query_paired = promise_me(lock, paired_queries)
        query_1 = promise_me(module_get_first, query_paired)
        query_2 = promise_me(module_get_second, query_paired)
        seeds_1 = promise_me(seeding, fm_pledge, query_1)
        seeds_2 = promise_me(seeding, fm_pledge, query_2)
        empty_1 = promise_me(collector, query_1, seeds_1)
        empty_2 = promise_me(collector, query_2, seeds_2)
        empty = promise_me(join_module, empty_1, empty_2)
        res_vec.append(promise_me(
            UnLock(parameter_manager, query_paired), empty))
    # compute all seeds
    res_vec.simultaneous_get(NUM_THREADS)
    del database

    # setup output dictionary structure
    read_background_data = []
    read_foreground_data = []
    seed_lines_forward_data = []
    seed_lines_reverse_complement_data = []
    sv_line_data = []
    sv_arrow_data = []
    sv_inversion_arrow_data = []

    section_list = []
    # get all sections that we have to render
    for nuc_seq, segments in collector.cpp_module.collection:
        for seeds in for_seed_sections(segments, fm_index, len(nuc_seq)):
            # compute very first position of section
            start_ref = min(
                x.start_ref if x.on_forward_strand else x.start_ref - x.size for x in seeds)
            section_list.append((start_ref, seeds, nuc_seq))
    # sort sections by start position on reference
    section_list.sort(key=lambda x: x[0])

    # total offset of graphic
    x_offset = section_list[0][0]

    # end_list is used to determine the y coordinate of the reads (checking for the overlaps...)
    end_list = []
    for start_ref, seeds, nuc_seq in section_list:
        # compute very last position of section
        end_ref = max(
            x.start_ref + x.size if x.on_forward_strand else x.start_ref for x in seeds)
        # compute the row id for the current sections so that there are no overlaps
        idx = len(end_list)
        for index, end_pos in enumerate(end_list):
            if end_pos + 3 < start_ref:
                idx = index
                break
        # if section would overlap ov every existant row, add new row
        if idx == len(end_list):
            end_list.append(0)
        # update the end list
        end_list[idx] = end_ref

        # add background panel for seed section
        read_background_data.append(
            [start_ref - x_offset, idx, end_ref - start_ref, 1, 1, str(nuc_seq.id)])

        # add foreground panels for seeds and seeds themselves
        seeds.sort(key=lambda x: x.start) # this is so that we have the seeds indexed according to query positions
        for seed in seeds:
            x = None
            if seed.on_forward_strand:
                x = seed.start_ref - x_offset
            else:
                x = seed.start_ref - x_offset - seed.size
            read_foreground_data.append([x, idx, seed.size, 1])

            if seed.on_forward_strand:
                seed_lines_forward_data.append([seed.start_ref - x_offset,
                                                seed.start / len(nuc_seq) + idx, seed.size, seed.size / len(nuc_seq)])
            else:
                seed_lines_reverse_complement_data.append([
                    seed.start_ref - x_offset,
                    seed.start / len(nuc_seq) + idx,
                    -seed.size,
                    seed.size / len(nuc_seq)
                ])

    conn = sqlite3.connect(db_name)
    cur = conn.cursor()
    # get all the sv jumps we have to render
    cur.execute(
        """ SELECT sv_line.start, sv_line.end, sv_jump.start, sv_jump.end, sv_jump.switch_strand
            FROM sv_jump
            JOIN sv_line ON sv_line.id = sv_jump.sv_line_id
    """)
    y = 0
    for a_start, a_end, b_start, b_end, switch_strand in cur.fetchall():
        print(a_start, a_end, b_start, b_end, switch_strand)
        ref_from = a_end
        ref_to = b_start
        if b_start < a_end:
            ref_from = a_end
            ref_to = b_start
        if switch_strand:
            sv_inversion_arrow_data.append(
                [ref_from - x_offset, y, ref_to - ref_from, 0])
        else:
            sv_arrow_data.append(
                [ref_from - x_offset, y, ref_to - ref_from, 0])
        y -= 1
    # show SV lines
    cur.execute(
        """ SELECT id, start, end
            FROM sv_line
    """)

    for idx, start_on_ref, end_on_ref in cur.fetchall():
        sv_line_data.append([start_on_ref - x_offset, 0,
                             end_on_ref - start_on_ref + 1, y])
    conn.close()

    # combine to single dictionary
    out_dict = {
        "x_offset": x_offset,
        "panels": [
            {
                "items": [
                    {
                        "type": "box-alpha",
                        "color": "#cccccc",
                        "line_color": None,
                        "line_width": 0,
                        "group": "read_background",
                        "data": read_background_data
                    },
                    {
                        "type": "box",
                        "color": "#595959",
                        "group": "read_foreground",
                        "data": read_foreground_data
                    },
                    {
                        "type": "line",
                        "color": "blue",
                        "group": "seed_lines_forward",
                        "data": seed_lines_forward_data
                    },
                    {
                        "type": "line",
                        "color": "orange",
                        "group": "seed_lines_reverse_complement",
                        "data": seed_lines_reverse_complement_data
                    }
                ],
                "h": 700
            },
            {
                "items": [
                    {
                        "type": "box",
                        "color": "purple",
                        "group": "sv_line",
                        "data": sv_line_data
                    },
                    {
                        "type": "arrow",
                        "color": "black",
                        "group": "sv_arrow",
                        "data": sv_arrow_data
                    },
                    {
                        "type": "arrow",
                        "color": "green",
                        "group": "sv_arrow",
                        "data": sv_inversion_arrow_data
                    }
                ],
                "h": 300
            }
        ]
    }
    return out_dict


if __name__ == "__main__":
    out_dict = create_json_from_db("/MAdata/databases/sv_simulated",
                                   "/MAdata/genome/human/GRCh38.p12/ma/genome")
    with open("/MAdata/tmp/sv_diagramm.json", "w") as json_out:
        json.dump(out_dict, json_out)
    #render_from_json("/MAdata/tmp/sv_diagramm.json", 7500000, 7550000)
