from MA import *

class SvLine:
    def __init__(self, start, end, desc):
        self.start = start
        self.end = end
        self.desc = desc
        self.id = None # will get updated once the SvLine is inserted into the Database

class SV_DB:
    def __init__(self, db_name, open_create_hint="open"):
        self.db = libMA.SV_DB(db_name, open_create_hint)

    def get_transaction_context(self, sv_caller_name, sv_caller_desc, parameter_set_manager=ParameterSetManager()):
        class SvInserter(libMA.SvInserter):
            def __init__(self, inserter):
                self.inserter = inserter

            def insert_sv_line(self, sv_line):
                sv_line.id = self.inserter.insert_sv_line(sv_line.start, sv_line.end, sv_line.desc)

            def connect_sv_lines(self, sv_line_from, sv_line_to, jump, desc):
                self.inserter.connect_sv_lines(sv_line_from.id, sv_line_to.id, jump, desc)

        return SvInserter(libMA.SvInserter(self.db, sv_caller_name, sv_caller_desc))

    def get_soc_inserter_module(self, parameter_set_manager=ParameterSetManager()):
        return libMA.SoCDbWriter(parameter_set_manager, libMA.SoCInserter(self.db))

    def get_reads_pledge(self, parameter_set_manager=ParameterSetManager(), single=True, paired=True):
        if single and paired:
            return promise_me(libMA.AllNucSeqFromSql(parameter_set_manager, self.db))
        if single:
            return promise_me(NucSeqFromSql(parameter_set_manager, self.db))
        if paired:
            return promise_me(PairedNucSeqFromSql(parameter_set_manager, self.db))

    def clear_soc_table(self):
        self.db.clear_soc_table()

    def has_socs(self):
        return self.db.has_socs()
