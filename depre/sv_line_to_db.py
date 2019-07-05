from sql_interface import *

##
# @brief overlaps sv lines and filters them
#


class SvLineToDb(VolatileModule):
    def __init__(self, sv_db, sv_inserter):
        VolatileModule.__init__(self)  # this line is crucial DON'T DELETE ME
        self.sv_db = sv_db
        self.sv_inserter = sv_inserter

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

class DummySvLineToDb(VolatileModule):
    def __init__(self, sv_db, sv_inserter):
        VolatileModule.__init__(self)  # this line is crucial DON'T DELETE ME
        self.sv_db = sv_db
        self.sv_inserter = sv_inserter

        self.insert_line = self.sv_inserter.insert_sv_line(0, 0, "dummy insert line")

        self.min_insert_size = 300

    ##
    # @brief
    # @details
    # Reimplemented from MA.aligner.Module.execute.
    def execute(self, input_vec):
        nuc_seq_id, soc, nuc_seq_len = input_vec

        sweep_list = []
        for seed in soc:
            sweep_list.append( (seed.start, True) )
            sweep_list.append( (seed.start + seed.size, False) )
        sweep_list.append( (nuc_seq_len, True) ) # so that we record a potential interval at the end of the read...

        assert False < True # necessary for tuple ordering
        sweep_list.sort()
        
        count_open_intervals = 0
        last_end_pos = 0 # so that we record a potential interval at the start of the read...
        
        for curr_pos, is_start in sweep_list:
            if is_start:
                if count_open_intervals == 0:
                    if curr_pos - last_end_pos > self.min_insert_size:
                        print("dummy insert connection:", nuc_seq_id, "from-to:", 
                                last_end_pos, curr_pos, "size", curr_pos - last_end_pos )
                        self.insert_line.insert_support( nuc_seq_id, last_end_pos, True )
                        self.insert_line.insert_support( nuc_seq_id, curr_pos, True )
                count_open_intervals += 1
            else:
                count_open_intervals -= 1
                if count_open_intervals == 0:
                    last_end_pos = curr_pos
        assert count_open_intervals == 1

        return None

class getSvLinesToDb:
    def __init__(self, db_name, parameter_manager):
        self.__db_nme = db_name
        self.__parameter_manager = parameter_manager
        self.sv_line_to_db = None
        self.dummy_sv_line_to_db = None
        self.__sv_db = None
        self.__sv_inserter = None

    def __enter__(self):
        self.__sv_db = SV_DB(self.__db_nme)
        print("clearing call tables")
        self.__sv_inserter = libMA.SvInserter(self.__sv_db.db, "caller_name", "desc")
        self.sv_line_to_db = SvLineToDb(self.__sv_db, self.__sv_inserter)
        #self.dummy_sv_line_to_db = DummySvLineToDb(self.__sv_db, self.__sv_inserter)

    def __exit__(self, exc_type, exc_val, exc_tb):
        #del self.dummy_sv_line_to_db
        del self.sv_line_to_db
        del self.__sv_inserter
        del self.__sv_db