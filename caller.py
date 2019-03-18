from soc_computer import *
from sv_line_overlapper import *
from sv_line_to_db import *

class Caller():
    def __init__(self, db_name, pack, fm_index, parameter_set_manager=ParameterSetManager()):
        self.db_name = db_name
        self.pack = pack
        self.fm_index = fm_index
        self.parameter_set_manager = parameter_set_manager
        self.parameter_set_manager.set_selected("PacBio")

        # compute SoC's if it hasn't been done yet
        compute_socs(self.sv_db, self.pack, self.fm_index, self.parameter_set_manager, force=False)
        self.sv_db = SV_DB(db_name)
        self.sv_db.db.clear_calls_table()

        add_sv_line_params(self.parameter_set_manager)
        add_sv_overlap_params(self.parameter_set_manager)

        self.sv_line_module = SvLineOverlapper(self.db_name,
                            self.pack, self.fm_index, self.parameter_set_manager)
        self.sv_line_saver = SvLineToDb(self.db_name, self.parameter_set_manager)

    def call(self):
        while not self.sv_line_module.is_finished():
            x = self.sv_line_module.execute(None)
            print(x.ref_pos_start, "-", x.ref_pos_end, ":", str(type(x)).split("'")[1])
            self.sv_line_saver.execute([x])


if __name__ == "__main__":
    fm_index = FMIndex()
    fm_index.load("/MAdata/genome/human/GRCh38.p12/ma/genome")
    pack = Pack()
    pack.load("/MAdata/genome/human/GRCh38.p12/ma/genome")

    caller = Caller("/MAdata/databases/sv_simulated", pack, fm_index)
    caller.call()