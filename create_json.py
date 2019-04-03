from MA import *
import json
import sqlite3
from render_json import render_from_json
NUM_THREADS = 32

def for_seed_sections(segments, fm_index, nuc_seq_len):
    seeds = [x for x in segments.extract_seeds(fm_index, 100, 0, nuc_seq_len, True)]
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

def create_json_from_db(db_name, pack_path, out_file_name):
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

    # setup output dictionary structure
    read_background_item = {
        "type": "box",
        "color": "#cccccc",
        "group": "read_background",
        "y": [],
        "x": [],
        "w": [],
        "h": 1
    }
    read_foreground_item = {
        "type": "box",
        "color": "#595959",
        "group": "read_foreground",
        "y": [],
        "x": [],
        "w": [],
        "h": 1
    }
    seed_lines_forward_item = {
        "type": "line",
        "color": "blue",
        "group": "seed_lines_forward",
        "y": [],
        "x": [],
        "w": [],
        "h": []
    }
    seed_lines_reverse_complement_item = {
        "type": "line",
        "color": "orange",
        "group": "seed_lines_reverse_complement",
        "y": [],
        "x": [],
        "w": [],
        "h": []
    }
    sv_line_item = {
        "type": "box",
        "color": "purple",
        "group": "sv_line",
        "y": 0,
        "x": [],
        "w": [],
        "h": 0
    }
    sv_arrow_item = {
        "type": "arrow",
        "color": "black",
        "group": "sv_arrow",
        "y": [],
        "x": [],
        "w": [],
        "h": 0
    }
    sv_inversion_arrow_item = {
        "type": "arrow",
        "color": "green",
        "group": "sv_arrow",
        "y": [],
        "x": [],
        "w": [],
        "h": 0
    }

    section_list = []
    # get all sections that we have to render
    for nuc_seq, segments in collector.cpp_module.collection:
        for seeds in for_seed_sections(segments, fm_index, len(nuc_seq)):
            # compute very first position of section
            start_ref = min(x.start_ref if x.on_forward_strand else x.start_ref - x.size for x in seeds)
            section_list.append( (start_ref, seeds, nuc_seq) )
    # sort sections by start position on reference
    section_list.sort(key=lambda x: x[0])

    # total offset of graphic
    x_offset = section_list[0][0]

    # end_list is used to determine the y coordinate of the reads (checking for the overlaps...)
    end_list = []
    for start_ref, seeds, nuc_seq in section_list:
        # compute very last position of section
        end_ref = max(x.start_ref + x.size if x.on_forward_strand else x.start_ref for x in seeds)
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
        read_background_item["x"].append(start_ref - x_offset)
        read_background_item["y"].append(idx)
        read_background_item["w"].append(end_ref - start_ref)

        # add foreground panels for seeds and seeds themselves
        for seed in seeds:
            if seed.on_forward_strand:
                read_foreground_item["x"].append(seed.start_ref - x_offset)
            else:
                read_foreground_item["x"].append(seed.start_ref - x_offset - seed.size)
            read_foreground_item["w"].append(seed.size)
            read_foreground_item["y"].append(idx)

            if seed.on_forward_strand:
                seed_lines_forward_item["x"].append(seed.start_ref - x_offset)
                seed_lines_forward_item["y"].append(seed.start / len(nuc_seq) + idx)
                seed_lines_forward_item["w"].append(seed.size)
                seed_lines_forward_item["h"].append(seed.size / len(nuc_seq))
            else:
                seed_lines_reverse_complement_item["x"].append(seed.start_ref - x_offset)
                seed_lines_reverse_complement_item["y"].append(seed.start / len(nuc_seq) + idx)
                seed_lines_reverse_complement_item["w"].append(-seed.size)
                seed_lines_reverse_complement_item["h"].append(seed.size / len(nuc_seq))

    conn = sqlite3.connect(db_name)
    cur = conn.cursor()
    # show SV lines
    cur.execute(
        """ SELECT id, start, end
            FROM sv_line
    """)

    for idx, start_on_ref, end_on_ref in cur.fetchall():
        sv_line_item['x'].append(start_on_ref - x_offset)
        sv_line_item['w'].append(end_on_ref - start_on_ref + 1)
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
            sv_inversion_arrow_item["x"].append(ref_from - x_offset)
            sv_inversion_arrow_item["w"].append(ref_to - ref_from)
            sv_inversion_arrow_item["y"].append(y)
        else:
            sv_arrow_item["x"].append(ref_from - x_offset)
            sv_arrow_item["w"].append(ref_to - ref_from)
            sv_arrow_item["y"].append(y)
        y -= 1
    sv_line_item['h'] = y
    conn.close()

    # combine to single dictionary
    out_dict = {
        "x_offset": x_offset,
        "panels": [
            {
                "items": [
                    read_background_item,
                    read_foreground_item,
                    seed_lines_forward_item,
                    seed_lines_reverse_complement_item
                ],
                "h": 700
            },
            {
                "items": [
                    sv_line_item,
                    sv_arrow_item,
                    sv_inversion_arrow_item
                ],
                "h": 300
            }
        ]
    }
    with open(out_file_name, "w") as json_out:
        json.dump(out_dict, json_out)



if __name__ == "__main__":
    create_json_from_db("/MAdata/databases/sv_simulated",
                        "/MAdata/genome/human/GRCh38.p12/ma/genome",
                        "/MAdata/tmp/sv_diagramm.json")
    render_from_json("/MAdata/tmp/sv_diagramm.json", 7500000, 7550000)
