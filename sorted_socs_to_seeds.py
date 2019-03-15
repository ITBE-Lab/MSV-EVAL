from sorted_soc_from_sql import *

##
# @brief returns sorted seeds given sorted soc's
#
class SortedSoCtoSortedSeeds(VolatileModule):
    def __init__(self, db_name, pack, fm_index, parameter_manager):
        VolatileModule.__init__(self)  # this line is crucial DON'T DELETE ME
        self.heap = heap.Heap(
            lt_comp=lambda x, y: x[0][x[2]].start_ref < y[0][y[2]].start_ref)

        self.helper = SortedSocFromSQl(
            db_name, pack, fm_index, parameter_manager)
        self.next_soc = None
        if not self.helper.is_finished():
            self.next_soc = self.helper.execute(None)
        self.re_fill_heap()

    def re_fill_heap(self):
        while not self.next_soc is None and \
                (self.heap.empty() or
                 self.next_soc[0][0].start_ref <= self.heap.peek()[0][self.heap.peek()[2]].start_ref):

            # push the actual SoC and its id
            self.heap.push((self.next_soc[0], self.next_soc[1], 0, self.next_soc[2]))
            self.next_soc = None
            if not self.helper.is_finished():
                self.next_soc = self.helper.execute(None)
        if self.heap.empty():
            self.set_finished()

    ##
    # @brief
    # @details
    # Reimplemented from MA.aligner.Module.execute.
    def execute(self, input_vec):
        soc, soc_id, seed_index, nuc_seq_len = self.heap.pop()
        if seed_index + 1 < len(soc):
            self.heap.push((soc, soc_id, seed_index+1, nuc_seq_len))
        self.re_fill_heap()
        return soc, seed_index, soc_id, nuc_seq_len


def test_SortedSoCtoSortedSeeds(fm_index, pack):
    print("testing SortedSoCtoSortedSeeds...")
    x = SortedSoCtoSortedSeeds("/MAdata/databases/sv_simulated",
                         pack, fm_index, ParameterSetManager())
    last = 0
    heap_size = 0
    num_steps = 0
    while not x.is_finished():
        soc, seed_index, _, _ = x.execute(None)
        this = soc[seed_index].start_ref
        num_steps += 1
        heap_size += len(x.heap)
        assert this >= last
        last = this
    print("success! Average heap size:", heap_size/num_steps)


if __name__ == "__main__":
    print("loading index...")
    fm_index = FMIndex()
    fm_index.load("/MAdata/genome/human/GRCh38.p12/ma/genome")
    pack = Pack()
    pack.load("/MAdata/genome/human/GRCh38.p12/ma/genome")
    print("done")

    test_SortedSocFromSQl(fm_index, pack)
    test_SortedSoCtoSortedSeeds(fm_index, pack)

    print("done")
