from sql_interface import *

##
# @brief overlaps sv lines and filters them
#


class SvLineToDb(VolatileModule):
    def __init__(self, db_name, parameter_manager):
        VolatileModule.__init__(self)  # this line is crucial DON'T DELETE ME

        self.sv_db = SV_DB(db_name)
        print("clearing call tables")
        self.sv_db.db.clear_calls_table()
        self.sv_inserter = libMA.SvInserter(self.sv_db.db, "caller_name", "desc")

    def close(self):
        del self.sv_inserter

    ##
    # @brief
    # @details
    # Reimplemented from MA.aligner.Module.execute.

    def execute(self, input_vec):
        sv_line = input_vec[0]

        line_context = self.sv_inserter.insert_sv_line(
            sv_line.ref_pos_start, sv_line.ref_pos_end, str(type(sv_line)).split("'")[1])

        for read_id, q_pos, forw_strand in sv_line.supporting:
            line_context.insert_support(read_id, q_pos, forw_strand)

        return None
