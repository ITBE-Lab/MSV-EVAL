from sorted_socs_to_seeds import *
from sv_line import *
from sorted_socs_to_seeds import *

##
# @brief returns sorted SV lines given sorted seeds
#
COMPUTE_OVERLAP_LINES = False


class SeedsToSVLines(VolatileModule):
    def __init__(self, db_name, pack, fm_index, parameter_manager):
        VolatileModule.__init__(self)  # this line is crucial DON'T DELETE ME
        self.max_sv_line_size = parameter_manager.get_selected().by_name(
            "max sv line size").get()
        self.min_sv_line_size = parameter_manager.get_selected().by_name(
            "min sv line size").get()
        self.sv_line_size_slope = parameter_manager.get_selected().by_name(
            "sv line size slope").get()
        self.sv_line_fuzziness = parameter_manager.get_selected().by_name(
            "sv line fuzziness").get()

        self.heap = heap.Heap(
            lt_comp=lambda x, y: x.ref_pos_start < y.ref_pos_start)

        self.helper = SortedSoCtoSortedSeeds(
            db_name, pack, fm_index, parameter_manager)

        # soc, seed_index, soc_id, nuc_seq_len = self.helper.execute(None)
        self.next_seed = None
        if not self.helper.is_finished():
            self.next_seed = self.helper.execute(None)
        self.re_fill_heap()
        self.unpacked_size = pack.unpacked_size()

    def dist_interval(self, a_s, a_e, b_s, b_e):
        if b_s <= a_e and a_s <= b_e:
            return 0  # intervals are overlapping
        return min(abs(a_e - b_s), abs(a_e - b_s))

    def sv_line_size(self, seed, next_seed):
        q_dist = self.dist_interval(seed.start, seed.start + seed.size,
                                    next_seed.start, next_seed.start + next_seed.size)
        r_dist = self.dist_interval(seed.start_ref, seed.start_ref + seed.size,
                                    next_seed.start_ref, next_seed.start_ref + next_seed.size)

        dist = max(0, max(q_dist, r_dist) - self.min_sv_line_size) * \
            self.sv_line_size_slope + self.min_sv_line_size
        return min(self.max_sv_line_size, int(dist))

    def to_sv_lines(self, soc, index, soc_id, query_len):
        seed = soc[index]
        if seed.size > self.sv_line_fuzziness * 2:
            yield SeedSvLine(soc_id, seed.start_ref + self.sv_line_fuzziness, seed.start_ref + seed.size - self.sv_line_fuzziness, None, seed.on_forward_strand, query_len)

        if index > 0:
            prev = soc[index-1]
            yield GapEndSvLine(soc_id, seed.start_ref - self.sv_line_size(seed, prev),
                               seed.start_ref + self.sv_line_fuzziness, seed.start, seed.on_forward_strand, query_len)

            # overlap line with previous seed...
            if seed.start_ref < prev.start_ref + prev.size and seed.on_forward_strand == prev.on_forward_strand \
                    and COMPUTE_OVERLAP_LINES:
                overlap_size = self.max_sv_line_size  # @todo
                yield OverlapSvLine(soc_id, prev.start_ref + prev.size, prev.start_ref + prev.size + overlap_size,
                                    prev.start + prev.size, seed.on_forward_strand, query_len)
        elif seed.start > 0:  # if the seed does not reach the beginning of the read
            prev_ = libMA.Seed(0, 0, seed.start_ref, seed.on_forward_strand)
            yield GapEndSvLine(soc_id, seed.start_ref - self.sv_line_size(seed, prev_),
                               seed.start_ref + self.sv_line_fuzziness, seed.start, seed.on_forward_strand, query_len)

        if index < len(soc) - 1:
            next_ = soc[index+1]
            yield GapStartSvLine(soc_id, seed.start_ref + seed.size - self.sv_line_fuzziness,
                                 seed.start_ref + seed.size + self.sv_line_size(seed, next_), seed.start + seed.size,
                                 seed.on_forward_strand, query_len)

            # overlap line with next seed...
            if next_.start_ref < seed.start_ref + seed.size and next_.on_forward_strand == seed.on_forward_strand \
                    and COMPUTE_OVERLAP_LINES:
                overlap_size = self.max_sv_line_size  # @todo
                yield OverlapSvLine(soc_id, next_.start_ref - overlap_size, next_.start_ref,
                                    next_.start, seed.on_forward_strand, query_len)

        elif seed.start + seed.size < query_len:  # if the seed does not reach the end of the read
            next_ = libMA.Seed(query_len, 0, seed.start_ref + seed.size, seed.on_forward_strand)
            yield GapStartSvLine(soc_id, seed.start_ref + seed.size - self.sv_line_fuzziness,
                                 seed.start_ref + seed.size + self.sv_line_size(seed, next_), seed.start + seed.size,
                                 seed.on_forward_strand, query_len)

    def re_fill_heap(self):
        while not self.next_seed is None and \
                (self.heap.empty() or
                 self.next_seed[0][self.next_seed[1]].start_ref -
                    self.max_sv_line_size <= self.heap.peek().ref_pos_start):

            # push the actual SoC and its id
            for sv_line in self.to_sv_lines(*self.next_seed):
                self.heap.push(sv_line)
            self.next_seed = None
            if not self.helper.is_finished():
                self.next_seed = self.helper.execute(None)
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


def add_sv_line_params(parameter_manager):
    parameter_manager.get_selected().register_parameter(libMA.AlignerParameterInt(
        "max sv line size", "maximal size of sv line", 6, "Structural Variants Caller", 0))  # 25
    parameter_manager.get_selected().register_parameter(libMA.AlignerParameterInt(
        "min sv line size", "minimal size of sv line", 6, "Structural Variants Caller", 0))
    parameter_manager.get_selected().register_parameter(libMA.AlignerParameterInt(
        "sv line fuzziness", "fuzziness of sv line", 6, "Structural Variants Caller", 3))
    parameter_manager.get_selected().register_parameter(libMA.AlignerParameterDouble(
        "sv line size slope", "...", 6, "Structural Variants Caller", 0.05))


def test_SeedsToSVLines(fm_index, pack):
    print("testing SeedsToSVLines...")
    parameter_manager = ParameterSetManager()
    add_sv_line_params(parameter_manager)

    x = SeedsToSVLines("/MAdata/databases/sv_simulated",
                       pack, fm_index, parameter_manager)

    last = 0
    heap_size = 0
    num_steps = 0
    while not x.is_finished():
        sv_line = x.execute(None)
        this = sv_line.ref_pos_start
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

    #test_SortedSocFromSQl(fm_index, pack)
    #test_SortedSoCtoSortedSeeds(fm_index, pack)
    test_SeedsToSVLines(fm_index, pack)

    print("done")
