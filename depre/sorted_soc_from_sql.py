import sqlite3
from MA import *
import heap

##
# @brief extracts fm-index intervals from sql
# @details
class SegmentVectorsFromSQL(VolatileModule):
    def __init__(self, db_name, fm_index, parameter_set_manager):
        VolatileModule.__init__(self)  # this line is crucial DON'T DELETE ME
        self.conn = sqlite3.connect(db_name)
        self.cur = self.conn.cursor()
        self.cur.execute("""
            SELECT read_table.sequence, read_table.id
            FROM read_table
        """)
        self.next = self.cur.fetchone()
        if self.next is None:
            print("No Data")
            self.set_finished()

        self.fm_index = fm_index
        self.seeding_module = BinarySeeding(parameter_set_manager)

    ##
    # @brief
    # @details
    # Reimplemented from MA.aligner.Module.execute.
    def execute(self, input_vec):
        nuc_seq_blob, nuc_seq_id = self.next
        nuc_seq = nuc_seq_from_bytes(nuc_seq_blob)

        seeds = self.seeding_module.execute(self.fm_index, nuc_seq)

        self.next = self.cur.fetchone()
        if self.next is None:
            self.set_finished()

        return nuc_seq_id, seeds, len(nuc_seq)


##
# @brief extracts SoC's from database
# @details
# returns them semi sorted -> by the order given by soc_start in the database
class AllSoCsOfRead(VolatileModule):
    def __init__(self, db_name, pack, fm_index, parameter_set_manager):
        VolatileModule.__init__(self)  # this line is crucial DON'T DELETE ME
        self.max_padding = parameter_set_manager.get_selected().by_name(
            "re seeding padding").get()
        self.num_socs = parameter_set_manager.get_selected().by_name(
            "Number of SoC's for SV calling").get()
        print("Number of SoC's for SV calling:", self.num_socs)
        self.conn = sqlite3.connect(db_name)
        self.cur = self.conn.cursor()
        self.cur.execute("""
            SELECT read_table.sequence, read_table.id
            FROM read_table
        """)
        self.next = self.cur.fetchone()
        if self.next is None:
            print("No Data")
            self.set_finished()

        self.pack = pack
        self.fm_index = fm_index
        self.seeding_module = BinarySeeding(parameter_set_manager)
        self.soc_module = StripOfConsideration(parameter_set_manager)
        self.fill_module = FillSeedSet(parameter_set_manager)

    ##
    # @brief
    # @details
    # Reimplemented from MA.aligner.Module.execute.
    def execute(self, input_vec):
        nuc_seq_blob, nuc_seq_id = self.next
        nuc_seq = nuc_seq_from_bytes(nuc_seq_blob)

        seeds = self.seeding_module.execute(self.fm_index, nuc_seq)
        socs = self.soc_module.execute(
            seeds, nuc_seq, self.pack, self.fm_index)

        all_seeds = libMA.Seeds()
        for _ in range(self.num_socs):
            if socs.empty():
                break
            seeds = socs.pop()

            filled_seeds = self.fill_module.execute(
                seeds, nuc_seq, self.fm_index, self.pack)
            all_seeds.extend(filled_seeds)

        self.next = self.cur.fetchone()
        if self.next is None:
            self.set_finished()

        return nuc_seq_id, all_seeds, len(nuc_seq)



##
# @brief extracts SoC's from database
# @details
# returns them semi sorted -> by the order given by soc_start in the database
class SoCSortedSocFromSQl(VolatileModule):
    def __init__(self, db_name, pack, fm_index, parameter_set_manager):
        VolatileModule.__init__(self)  # this line is crucial DON'T DELETE ME
        self.max_padding = parameter_set_manager.get_selected().by_name(
            "re seeding padding").get()
        self.conn = sqlite3.connect(db_name)
        self.cur = self.conn.cursor()
        self.cur.execute("""
            SELECT read_table.sequence, read_table.id, soc_table.soc_index, soc_table.soc_start
            FROM read_table
            JOIN soc_table ON soc_table.read_id == read_table.id
            ORDER BY soc_table.soc_start
        """)
        self.next = self.cur.fetchone()
        if self.next is None:
            print("No Data")
            self.set_finished()

        self.pack = pack
        self.fm_index = fm_index
        self.seeding_module = BinarySeeding(parameter_set_manager)
        self.soc_module = StripOfConsideration(parameter_set_manager)
        self.fill_module = FillSeedSet(parameter_set_manager)

    ##
    # @brief
    # @details
    # Reimplemented from MA.aligner.Module.execute.
    def execute(self, input_vec):
        nuc_seq_blob, nuc_seq_id, soc_index, soc_start = self.next
        nuc_seq = nuc_seq_from_bytes(nuc_seq_blob)

        seeds = self.seeding_module.execute(self.fm_index, nuc_seq)
        socs = self.soc_module.execute(
            seeds, nuc_seq, self.pack, self.fm_index)

        assert len(socs) > soc_index
        for _ in range(soc_index):
            socs.pop()
        seeds = socs.pop()

        filled_seeds = self.fill_module.execute(
            seeds, nuc_seq, self.fm_index, self.pack)

        # sort filled seeds by reference position and assert that it worked
        filled_seeds.sort_by_ref_pos()
        for a, b in zip(filled_seeds[:-1], filled_seeds[1:]):
            assert a.start_ref <= b.start_ref
        if soc_start - self.max_padding > filled_seeds[0].start_ref:
            print("soc_start:", soc_start, "actual_start",
                  filled_seeds[0].start_ref, "diff:", soc_start -
                  filled_seeds[0].start_ref,
                  "max_allowed_diff:", self.max_padding)
            exit(0)

        self.next = self.cur.fetchone()
        if self.next is None:
            self.set_finished()

        return min(soc_start, filled_seeds[0].start_ref), (nuc_seq_id, soc_index), filled_seeds, len(nuc_seq)


##
# @brief orders SoC's extracted by SoCSortedSocFromSQl
# @details
# orders the soc's by the first seeds reference position
class SortedSocFromSQl(VolatileModule):
    def __init__(self, db_name, pack, fm_index, parameter_manager):
        VolatileModule.__init__(self)  # this line is crucial DON'T DELETE ME
        self.max_padding = parameter_manager.get_selected().by_name(
            "re seeding padding").get()

        self.helper = SoCSortedSocFromSQl(
            db_name, pack, fm_index, parameter_manager)
        self.heap = heap.Heap(
            lt_comp=lambda x, y: x[0][0].start_ref < y[0][0].start_ref)

        self.next_soc = None
        if not self.helper.is_finished():
            self.next_soc = self.helper.execute(None)

        self.re_fill_heap()

    def re_fill_heap(self):
        while not self.next_soc is None and \
            (self.heap.empty() or
             self.next_soc[0] - self.max_padding <= self.heap.peek()[0][0].start_ref):

            # push the actual SoC and its id
            self.heap.push((self.next_soc[2], self.next_soc[1], self.next_soc[3]))
            if self.helper.is_finished():
                self.next_soc = None
            else:
                self.next_soc = self.helper.execute(None)
        if self.heap.empty():
            self.set_finished()

    ##
    # @brief
    # @details
    # Reimplemented from MA.aligner.Module.execute.
    def execute(self, input_vec):
        soc, soc_id, nuc_seq_len = self.heap.pop()
        self.re_fill_heap()
        return soc, soc_id, nuc_seq_len

def test_SortedSocFromSQl(fm_index, pack):
    print("testing SortedSocFromSQl...")
    x = SortedSocFromSQl("/MAdata/databases/sv_simulated",
                         pack, fm_index, ParameterSetManager())
    last = 0
    heap_size = 0
    num_steps = 0
    while not x.is_finished():
        this = x.execute(None)[0][0].start_ref
        num_steps += 1
        heap_size += len(x.heap)
        assert this >= last
        last = this
    print("success! Average heap size:", heap_size/num_steps)