from MA import *
from compute_sv_jumps import *
from sweep_sv_jumps import *
from render_json import render_from_dict
from create_json import create_json_from_db
import sqlite3

class Caller():
    def __init__(self, db_name, pack, fm_index, parameter_set_manager=ParameterSetManager()):
        self.db_name = db_name
        self.conn = sqlite3.connect(db_name)
        self.pack = pack
        self.fm_index = fm_index
        self.parameter_set_manager = parameter_set_manager
        self.parameter_set_manager.set_selected("PacBio")
        self.parameter_set_manager.get_selected().register_parameter(
            # name, description, group_id, group_name, value
            libMA.AlignerParameterInt("num destinations", "desc", 6, "Structural Variants Caller", 5))
        self.parameter_set_manager.get_selected().register_parameter(
            # name, description, group_id, group_name, value
            libMA.AlignerParameterInt("fuzziness", "desc", 6, "Structural Variants Caller", 3))
        self.parameter_set_manager.get_selected().register_parameter(
            # name, description, group_id, group_name, value
            libMA.AlignerParameterInt("min coverage", "desc", 6, "Structural Variants Caller", 7))

    def close(self):
        self.conn.commit()
        self.conn.close()

    def get_jumps(self):
        return compute_sv_jumps(self.parameter_set_manager, self.conn, self.fm_index)

    def call(self):
        print("calling sv lines...")
        sv_jumps = self.get_jumps()
        print("sweeping sv lines...")
        accepted_sv_jumps = sweep_sv_jumps(self.parameter_set_manager, self.conn, sv_jumps)
        print("done")


if __name__ == "__main__":
    fm_index = FMIndex()
    fm_index.load("/MAdata/genome/human/GRCh38.p12/ma/genome")
    pack = Pack()
    pack.load("/MAdata/genome/human/GRCh38.p12/ma/genome")

    caller = Caller("/MAdata/databases/sv_simulated", pack, fm_index)

    out_dict = sv_jumps_to_dict(caller.get_jumps())
    render_from_dict(out_dict, 7500000, 7550000, True)
    caller.close()
    exit()

    caller.call()
    caller.close()
    
    # display the result
    out_dict = create_json_from_db("/MAdata/databases/sv_simulated",
                    "/MAdata/genome/human/GRCh38.p12/ma/genome")
    render_from_dict(out_dict, 7500000, 7550000)
