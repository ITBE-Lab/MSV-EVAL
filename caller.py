from soc_computer import *
from sv_line_overlapper import *
from sv_line_to_db import *
from extract_connected_sv_lines import *
from sql_interface import *

class Caller():
    def __init__(self, db_name, pack, fm_index, parameter_set_manager=ParameterSetManager()):
        self.db_name = db_name
        self.pack = pack
        self.fm_index = fm_index
        self.parameter_set_manager = parameter_set_manager
        self.parameter_set_manager.set_selected("PacBio")

        # compute SoC's if it hasn't been done yet
        compute_socs(SV_DB(db_name), self.pack, self.fm_index, self.parameter_set_manager, force=False)

        add_sv_line_params(self.parameter_set_manager)
        add_sv_overlap_params(self.parameter_set_manager)

        self.sv_line_module = SvLineOverlapper(self.db_name,
                            self.pack, self.fm_index, self.parameter_set_manager)
        self.all_soc_from_sql = AllSoCsOfRead(self.db_name,
                            self.pack, self.fm_index, self.parameter_set_manager)
        self.sv_line_saver = getSvLinesToDb(self.db_name, self.parameter_set_manager)

        self.supported_extractor = ExtractSupportedSvLines(self.db_name, self.parameter_set_manager)
        self.line_connector = ConnectSvLines(self.db_name, self.parameter_set_manager)

    def call(self):
        print("calling sv lines...")
        i = 1
        with self.sv_line_saver:
            while not self.sv_line_module.is_finished():
                x = self.sv_line_module.execute(None)
                print(i, ";", x.ref_pos_start, "-", x.ref_pos_end, ":", str(type(x)).split("'")[1])
                i += 1
                self.sv_line_saver.sv_line_to_db.execute([x])
            #while not self.all_soc_from_sql.is_finished():
            #    x = self.all_soc_from_sql.execute(None)
            #    self.sv_line_saver.dummy_sv_line_to_db.execute(x)
        del self.sv_line_module
        del self.all_soc_from_sql
        print("done")

        print("connecting sv lines..")
        self.supported_extractor.load_ids()
        while not self.supported_extractor.is_finished():
            x = self.supported_extractor.execute(None)
            self.line_connector.execute([x])
        self.line_connector.close()
        print("done")

        print("filtering sv lines & connectors..")
        filter_sv_line_connectors(self.db_name)
        filter_sv_lines(self.db_name)
        print("done")


if __name__ == "__main__":
    fm_index = FMIndex()
    fm_index.load("/MAdata/genome/human/GRCh38.p12/ma/genome")
    pack = Pack()
    pack.load("/MAdata/genome/human/GRCh38.p12/ma/genome")

    caller = Caller("/MAdata/databases/sv_simulated", pack, fm_index)
    caller.call()