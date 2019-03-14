from soc_computer import *

class Caller():
    def __init__(self, db_name, pack, fm_index, parameter_set_manager=ParameterSetManager()):
        self.db_make = db_name
        self.sv_db = SV_DB(db_name)
        self.pack = pack
        self.fm_index = fm_index
        self.parameter_set_manager = parameter_set_manager

        # compute SoC's if it hasn't been done yet
        compute_socs(sv_db, pack, fm_index, parameter_set_manager, force=False)

    def call(self):
        pass


if __name__ == "__main__":
    fm_index = FMIndex()
    fm_index.load("/MAdata/genome/human/GRCh38.p12/ma/genome")
    pack = Pack()
    pack.load("/MAdata/genome/human/GRCh38.p12/ma/genome")

    caller = Caller("/MAdata/databases/sv_simulated", pack, fm_index)
    caller.call()