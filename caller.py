from MA import *
from compute_sv_jumps import *
from sweep_sv_jumps import *
from render_json import render_from_dict
from create_json import create_json_from_db
import sqlite3
import os
import json

class Caller():
    def __init__(self, db_name, pack, fm_index, parameter_set_manager=ParameterSetManager()):
        self.db_name = db_name
        self.pack = pack
        self.fm_index = fm_index
        self.parameter_set_manager = parameter_set_manager
        self.parameter_set_manager.set_selected("PacBio")
        self.parameter_set_manager.get_selected().register_parameter(
            # name, description, group_id, group_name, value
            libMA.AlignerParameterInt("num destinations", "desc", 6, "Structural Variants Caller", 2))
        self.parameter_set_manager.get_selected().register_parameter(
            # name, description, group_id, group_name, value
            libMA.AlignerParameterInt("fuzziness", "desc", 6, "Structural Variants Caller", 3))
        self.parameter_set_manager.get_selected().register_parameter(
            # name, description, group_id, group_name, value
            libMA.AlignerParameterInt("min coverage", "desc", 6, "Structural Variants Caller", 7))

    def call(self):
        print("calling sv jumps...")
        sv_db = SV_DB(self.db_name, "open")
        
        sv_db.clear_calls_table_for_caller("MA-SV") # do clear old runs?

        sv_db.set_num_threads(self.parameter_set_manager.get_num_threads())
        run_id = compute_sv_jumps(self.parameter_set_manager, self.fm_index, self.pack, sv_db)
        print("called", sv_db.num_jumps(), "jumps.")
        print("clustering sv jumps...")
        sweep_sv_jumps(self.parameter_set_manager, sv_db, run_id,
                                           self.pack.unpacked_size_single_strand)
        print("done")
        return sv_db, run_id


if __name__ == "__main__":
    fm_index = FMIndex()
    fm_index.load("/MAdata/genome/human/GRCh38.p12/ma/genome")
    pack = Pack()
    pack.load("/MAdata/genome/human/GRCh38.p12/ma/genome")

    caller = Caller("/MAdata/databases/sv_simulated", pack, fm_index)
    sv_db, run_id = caller.call()
    
    out_dict = sv_jumps_to_dict(sv_db, run_id)
    #with open("/MAdata/tmp/sv_diagramm.json", "w") as json_out:
    #    json.dump(out_dict, json_out)
    render_from_dict(out_dict, 7500000, 7550000, True)
    
    # display the result
    out_dict = create_json_from_db(sv_db, "/MAdata/genome/human/GRCh38.p12/ma/genome")
    render_from_dict(out_dict, 7500000, 7550000)
