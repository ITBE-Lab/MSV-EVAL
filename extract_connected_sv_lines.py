import sqlite3
from libMA import *

class SvLineConn:
    def __init__(self, soc_line_id, other_lines):
        self.soc_line_id = soc_line_id
        self.by_line_id = {}
        for read_id, read_pos, sv_line_id in other_lines:
            if sv_line_id not in self.by_line_id:
                self.by_line_id[sv_line_id] = {}
            if read_id not in self.by_line_id[sv_line_id]:
                self.by_line_id[sv_line_id][read_id] = []
            self.by_line_id[sv_line_id][read_id].append(read_pos)


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
            SELECT read_id, read_pos, sv_line_id
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

    def close(self):
        self.conn.commit()
        self.conn.close()

    def is_better_distance_list(self, current, alternative):
        print(current, "vs", alternative)
        if len(current) == 0:
            return True
        if len(current) < 3 and len(alternative) >= 3:
            return True
        if len(alternative) < 3 and len(current) >= 3:
            return False
        min_curr = min(current)
        min_alt = min(alternative)
        return min_alt < min_curr

    ##
    # @brief
    # @details
    # Reimplemented from MA.aligner.Module.execute.
    def execute(self, input_vec):
        sv_line_conn = input_vec[0]

        id_to = None
        to_dists = []

        for potential_match in sv_line_conn.by_line_id.keys():
            if potential_match == sv_line_conn.soc_line_id:
                continue
            dists = []
            for read_id, pos_to_s in sv_line_conn.by_line_id[potential_match].items():
                pos_from_s = sv_line_conn.by_line_id[sv_line_conn.soc_line_id][read_id]

                min_dist = abs(pos_from_s[0] - pos_to_s[0])
                for pos_from in pos_from_s:
                    for pos_to in pos_to_s:
                        min_dist = min(abs(pos_from - pos_to), min_dist)
                dists.append(min_dist)
            if self.is_better_distance_list(to_dists, dists):
                print("Better")
                id_to = potential_match
                to_dists = dists
            else:
                print("Worse")

        print("connecting", sv_line_conn.soc_line_id, "to", id_to)
        if not id_to is None:
            self.cur.execute("""
                INSERT INTO sv_line_connector_table
                VALUES (NULL, ?, ?, ?, ?)
            """, (sv_line_conn.soc_line_id, id_to, True, ""))

        return None

def filter_sv_line_connectors(db_name):
    conn = sqlite3.connect(db_name)
    cur = conn.cursor()
    sql_passing = """
        SELECT DISTINCT A.id
        FROM sv_line_table
        INNER JOIN sv_line_connector_table A ON A.soc_line_from = sv_line_table.id
        INNER JOIN sv_line_connector_table B ON B.soc_line_from = A.soc_line_to AND B.soc_line_to = sv_line_table.id
    """
    cur.execute(sql_passing)
    sv_line_conn_ids = cur.fetchall()
    print("interlaced connectors PASSED:", sv_line_conn_ids)
    cur.executemany("""
        INSERT INTO sv_line_connector_filter_table
        VALUES (NULL, ?, 1, "interlaced connectors")
    """, sv_line_conn_ids)

    cur.execute("""
        SELECT DISTINCT id
        FROM sv_line_connector_table
        WHERE id NOT IN (
        """ + sql_passing + ")")
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

