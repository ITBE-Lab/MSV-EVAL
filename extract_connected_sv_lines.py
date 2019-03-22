import sqlite3
from libMA import *
import math

class SvLineConn:
    def __init__(self, soc_line_id, other_lines):
        self.soc_line_id = soc_line_id
        self.by_line_id = {}
        for read_id, read_pos, sv_line_id, forw_strand in other_lines:
            if sv_line_id not in self.by_line_id:
                self.by_line_id[sv_line_id] = {}
            if read_id not in self.by_line_id[sv_line_id]:
                self.by_line_id[sv_line_id][read_id] = []
            self.by_line_id[sv_line_id][read_id].append( (read_pos, forw_strand) )


##
# @brief extracts sv lines that are connected via the sv_line_support_table
#
class ExtractSupportedSvLines(VolatileModule):
    def __init__(self, db_name, parameter_manager):
        VolatileModule.__init__(self)  # this line is crucial DON'T DELETE ME

        self.conn = sqlite3.connect(db_name)
        self.cur = self.conn.cursor()
        self.sv_line_ids = []

    def load_ids(self):
        self.cur.execute("SELECT id FROM sv_line_table")
        self.sv_line_ids = self.cur.fetchall()

        if len(self.sv_line_ids) == 0:
            self.set_finished()
            self.conn.close()

    ##
    # @brief
    # @details
    # Reimplemented from MA.aligner.Module.execute.
    def execute(self, input_vec):
        self.cur.execute("""
            SELECT read_id, read_pos, sv_line_id, on_forward_strand
            FROM sv_line_support_table
            WHERE read_id IN (
                SELECT read_id
                FROM sv_line_support_table
                WHERE sv_line_id = ?
            )
        """, self.sv_line_ids[0])
        ret = SvLineConn(self.sv_line_ids[0][0], self.cur.fetchall())
        if len(self.sv_line_ids) == 1:
            self.set_finished()
            self.conn.close()

        del self.sv_line_ids[0]
        return ret


##
# @brief 
#
class ConnectSvLines(Module):
    def __init__(self, db_name, parameter_manager):
        Module.__init__(self)  # this line is crucial DON'T DELETE ME
        self.conn = sqlite3.connect(db_name)
        self.cur = self.conn.cursor()
        self.max_jum_dist = 30
        self.score_count = 0.5
        self.penalty_deviation = 0.05

    def close(self):
        self.conn.commit()
        self.conn.close()

    def score_for_distance_list(self, dist_list):
        score = 0
        for x in dist_list:
            score += 1 / math.log(x + 1.5)
        return score

    # dynamic programming to find the best score
    # @todo this can be optimized for sure...
    def insertion_score(self, dist_list):
        dist_list.sort()
        max_score = 0
        switch_strand = False
        for start_pos in range(0, len(dist_list)):
            for end_pos in range(start_pos, len(dist_list)):
                median_dist, _ = dist_list[int((start_pos + end_pos) / 2)]
                total_deviation = sum(abs(x - median_dist) for x, _ in dist_list[start_pos:end_pos])*self.penalty_deviation
                total_count = (end_pos - start_pos)*self.score_count
                if total_count - total_deviation > max_score:
                    max_score = total_count - total_deviation
                    switch_strand = sum( 1 if x else -1 for _, x in dist_list[start_pos:end_pos]) > 0
        return max_score, switch_strand


        #if len(current) == 0:
        #    return True
        #if len(current) < 3 and len(alternative) >= 3:
        #    return True
        #if len(alternative) < 3 and len(current) >= 3:
        #    return False
        #current.sort()
        #alternative.sort()
        #min_curr = current[int(len(current)*0.1)]
        #min_alt = alternative[int(len(alternative)*0.1)]
        #if abs(min_curr - min_alt) < 10: # if the distance is less than 10 nt...
        #    return sum(alternative)/len(alternative) < sum(current)/len(current)
        #return min_alt < min_curr

    ##
    # @brief
    # @details
    # Reimplemented from MA.aligner.Module.execute.
    def execute(self, input_vec):
        sv_line_conn = input_vec[0]

        id_to = None
        to_dists = []
        switch_strand = None
        best_score = 0

        print("connecting:", sv_line_conn.soc_line_id)

        for potential_match in sv_line_conn.by_line_id.keys():
            if potential_match == sv_line_conn.soc_line_id: # if we are comparing with ourself -> skip
                continue
            print("checking with", potential_match)
            dists = []
            switch_strands = []
            for read_id, pos_to_s in sv_line_conn.by_line_id[potential_match].items():
                pos_from_s = sv_line_conn.by_line_id[sv_line_conn.soc_line_id][read_id]

                min_dist = None
                curr_switch_strand = None
                for pos_from, forw_strand_from in pos_from_s:
                    for pos_to, forw_strand_to in pos_to_s:
                        dist = abs(pos_from - pos_to)
                        if min_dist is None or dist < min_dist or (
                                dist == min_dist and forw_strand_from == forw_strand_to):
                            min_dist = dist
                            curr_switch_strand = forw_strand_from != forw_strand_to
                if not min_dist is None:
                    dists.append(min_dist)
                    switch_strands.append(curr_switch_strand)
            curr_score = self.score_for_distance_list(dists)
            print(to_dists, "->", best_score, "vs", dists, "->", curr_score)
            if curr_score > best_score:
                print("Better")
                id_to = potential_match
                best_score = curr_score
                to_dists = dists
                switch = 0
                no_switch = 0
                for s in switch_strands:
                    if s:
                        switch += 1
                    else:
                        no_switch += 1
                switch_strand = switch * 4 > no_switch
                print("switch strands? ", switch_strands, "->", "yes" if switch_strand else "no")
            else:
                print("Worse")

        if best_score < 5:
            print("checking with insertion")
            insertion_dists = []
            for pos_s in sv_line_conn.by_line_id[sv_line_conn.soc_line_id].values():
                for x, (pos_from, from_forw) in enumerate(pos_s):
                    for pos_to, to_forw in pos_s[0:x]:
                        if pos_from == pos_to:
                            continue
                        insertion_dists.append( (abs(pos_from - pos_to), from_forw != to_forw) )
            insertion_score, switch_strand_insertion = self.insertion_score(insertion_dists)
            print(insertion_dists, "->", insertion_score)

            if insertion_score > best_score:
                print("Better")
                id_to = sv_line_conn.soc_line_id
                switch_strand = switch_strand_insertion
                best_score = insertion_score
            else:
                print("Worse")

        print("connecting", sv_line_conn.soc_line_id, "to", id_to, "switch strand:", switch_strand, "final score:", best_score)
        if not id_to is None:
            self.cur.execute("""
                INSERT INTO sv_line_connector_table
                VALUES (NULL, ?, ?, ?, ?, ?)
            """, (sv_line_conn.soc_line_id, id_to, True, switch_strand, ""))

        return None

def filter_sv_line_connectors(db_name):
    conn = sqlite3.connect(db_name)
    cur = conn.cursor()
    # get all sv line connectors, that have a mate forming a two line cycle.
    sql_passing = """
        SELECT DISTINCT A.id
        FROM sv_line_table
        INNER JOIN sv_line_connector_table A ON A.soc_line_from = sv_line_table.id
        INNER JOIN sv_line_connector_table B ON B.soc_line_from = A.soc_line_to AND B.soc_line_to = sv_line_table.id
        WHERE A.switch_strand == B.switch_strand
    """
    cur.execute(sql_passing)
    sv_line_conn_ids = cur.fetchall()
    # and also get those that connect to dummy sv lines
    sql_passing_2 = """
        SELECT DISTINCT id
        FROM sv_line_connector_table
        WHERE soc_line_to == 1
    """
    cur.execute(sql_passing_2)
    sv_line_conn_ids.extend(cur.fetchall())

    print("interlaced connectors PASSED:", sv_line_conn_ids)
    cur.executemany("""
        INSERT INTO sv_line_connector_filter_table
        VALUES (NULL, ?, 1, "interlaced connectors")
    """, sv_line_conn_ids)

    cur.execute("""
        SELECT DISTINCT id
        FROM sv_line_connector_table
        WHERE id NOT IN (""" + sql_passing + ") AND id NOT IN (" + sql_passing_2 + ")")
    sv_line_conn_failed_ids = cur.fetchall()
    print("interlaced connectors FAILED:", sv_line_conn_failed_ids)
    cur.executemany("""
        INSERT INTO sv_line_connector_filter_table
        VALUES (NULL, ?, 0, "interlaced connectors")
    """, sv_line_conn_failed_ids)

    conn.commit()
    conn.close()

def filter_sv_lines(db_name):
    conn = sqlite3.connect(db_name)
    cur = conn.cursor()
    # get all sv lines that are the origin of at least one connector that passed all filters
    # (since we filtered connectors before, this will mean, that there are interlaced connectors for this sv line)
    sql_passing = """
        SELECT DISTINCT id
        FROM sv_line_table
        WHERE id IN (
            SELECT soc_line_from 
            FROM sv_line_connector_table
            WHERE id NOT IN (
                SELECT sv_line_connector_id
                FROM sv_line_connector_filter_table
                WHERE sv_line_connector_filter_table.pass = 0
            )
        )
    """
    cur.execute(sql_passing)
    sv_line_ids = cur.fetchall()
    print("connected sv lines PASSED:", sv_line_ids)
    cur.executemany("""
        INSERT INTO sv_line_filter_table
        VALUES (NULL, ?, 1, "connected sv lines")
    """, sv_line_ids)

    cur.execute("""
        SELECT DISTINCT id
        FROM sv_line_table
        WHERE id NOT IN (
        """ + sql_passing + ")")
    sv_line_failed_ids = cur.fetchall()
    print("connected sv lines FAILED:", sv_line_failed_ids)
    cur.executemany("""
        INSERT INTO sv_line_filter_table
        VALUES (NULL, ?, 0, "connected sv lines")
    """, sv_line_failed_ids)

    conn.commit()
    conn.close()

