from sorted_socs_to_seeds import *
from sv_line import *
from sorted_socs_to_seeds import *

##
# @brief returns sorted SV lines given sorted seeds
#


class SeedsToSVLines(VolatileModule):
    def __init__(self, db_name, pack, fm_index, parameter_manager):
        VolatileModule.__init__(self)  # this line is crucial DON'T DELETE ME
        self.max_sv_line_size = parameter_manager.get_selected().by_name("max sv line size")
        self.min_sv_line_size = parameter_manager.get_selected().by_name("min sv line size")
        self.sv_line_size_slope = parameter_manager.get_selected().by_name("sv line size slope")
        self.sv_line_fuzziness = parameter_manager.get_selected().by_name("sv line fuzziness")

        self.heap = heap.Heap(
            lt_comp=lambda x, y: x.ref_pos_start < y.ref_pos_start)

        self.helper = SortedSoCtoSortedSeeds(
            db_name, pack, fm_index, parameter_manager)

        # soc, seed_index, soc_id, nuc_seq_len = self.helper.execute(None)
        self.next_seed = None
        if not self.helper.is_finished():
            self.next_seed = self.helper.execute(None)
        self.re_fill_heap()

    def dist_interval(a_s, a_e, b_s, b_e):
        if b_s <= a_e and a_s <= b_e:
            return 0  # intervals are overlapping
        return min(abs(a_e - b_s), abs(a_e - b_s))

    def sv_line_size(self, seed, next_seed):
        q_dist = dist_interval(seed.start, seed.start + seed.size,
                               next_seed.start, next_seed.start + next_seed.size)
        r_dist = dist_interval(seed.start_ref, seed.start_ref + seed.size,
                               next_seed.start_ref, next_seed.start_ref + next_seed.size)

        dist = max(0, max(q_dist, r_dist) - self.min_sv_line_size) * \
                   self.sv_line_size_slope + self.min_sv_line_size
        return min(self.max_sv_line_size, dist)

    def to_sv_lines(self, soc, index, soc_id, query_len):
        seed = soc[index]
        if seed.size > self.sv_line_fuzziness * 2:
            yield SeedSvLine(soc_id, seed.start_ref + self.sv_line_fuzziness, seed.start_ref + seed.size - self.sv_line_fuzziness)

        # @todo continue here

    def re_fill_heap(self):
        while not self.next_seed is None and \
                (self.heap.empty() or
                 self.next_seed[0][self.next_seed[1]].start_ref] - \
                    self.max_sv_line_size <= self.heap.peek().ref_pos_start):

            # push the actual SoC and its id
            for sv_line in self.to_sv_lines(*self.next_seed):
                self.heap.push(sv_line)
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
        sv_line = self.heap.pop()
        self.re_fill_heap()
        return sv_line


def test_SeedsToSVLines(fm_index, pack):
    print("testing SeedsToSVLines...")
    parameter_manager= ParameterSetManager()
    parameter_manager.get_selected().register_parameter(AlignerParameterInt(
        "max sv line size", "maximal size of sv line", 6, "Structural Variants Caller", 25))
    parameter_manager.get_selected().register_parameter(AlignerParameterInt(
        "min sv line size", "minimal size of sv line", 6, "Structural Variants Caller", 3))
    parameter_manager.get_selected().register_parameter(AlignerParameterInt(
        "sv line fuzziness", "fuzziness of sv line", 6, "Structural Variants Caller", 3))
    parameter_manager.get_selected().register_parameter(AlignerParameterDouble(
        "sv line size slope", "...", 6, "Structural Variants Caller", 0.5))
    x = SortedSoCtoSortedSeeds("/MAdata/databases/sv_simulated",
                               pack, fm_index, parameter_manager)
    last=0
    heap_size=0
    num_steps=0
    while not x.is_finished():
        pass
        #soc, seed_index, _, _=x.execute(None)
        #this=soc[seed_index].start_ref
        #num_steps += 1
        #heap_size += len(x.heap)
        #assert this >= last
        #last=this
    print("success! Average heap size:", heap_size/num_steps)


if __name__ == "__main__":
    print("loading index...")
    fm_index=FMIndex()
    fm_index.load("/MAdata/genome/human/GRCh38.p12/ma/genome")
    pack=Pack()
    pack.load("/MAdata/genome/human/GRCh38.p12/ma/genome")
    print("done")

    test_SortedSocFromSQl(fm_index, pack)
    test_SortedSoCtoSortedSeeds(fm_index, pack)
    test_SeedsToSVLines(fm_index, pack)

    print("done")
