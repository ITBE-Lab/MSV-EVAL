from seeds_to_sv_lines import *

class OpenSvLineHeaps:
    def __init__(self):
        def lt_comp(x, y): return x.ref_pos_end < y.ref_pos_end
        self.heaps = {
            SeedSvLine: Heap(lt_comp),
            GapStartSvLine: Heap(lt_comp),
            GapEndSvLine: Heap(lt_comp),
            OverlapSvLine: Heap(lt_comp)
        }

    def push(self, sv_line):
        self.heaps[type(sv_line)].push(sv_line)

    def first(self):
        return min(self.heaps.values(), key=lambda x: x.peek().ref_pos_end)

    def peek(self):
        return self.first().peek()

    def pop(self):
        return self.first().pop()

    def empty():
        for x in self.heaps.values():
            if not x.empty():
                return False
        return True

    def num_contradicting(self, ref_pos_start):
        contradicting = set()
        for sv_line in self.heaps[SeedSvLine]:
            if sv_line.ref_pos_start <= ref_pos_start:
                contradicting.add(sv_line.nuc_seq_id)
        return len(contradicting)


##
# @brief overlaps sv lines and filters them
#
class SvLineOverlapper(VolatileModule):
    def __init__(self, db_name, pack, fm_index, parameter_manager):
        VolatileModule.__init__(self)  # this line is crucial DON'T DELETE ME

        self.heap = OpenSvLineHeaps()

        self.helper = SeedsToSVLines(
            db_name, pack, fm_index, parameter_manager)

        self.next_sv_line = None
        if not self.helper.is_finished():
            self.next_sv_line = self.helper.execute(None)
        self.re_fill_heap()


    def re_fill_heap(self):
        while not self.next_sv_line is None and \
                (self.heap.empty() or
                 self.next_sv_line.ref_pos_start <= self.heap.peek().ref_pos_end):

            # push the sv line
            self.heap.push(self.next_sv_line)
            self.next_sv_line = None
            if not self.helper.is_finished():
                self.next_sv_line = self.helper.execute(None)
        if self.heap.empty():
            self.set_finished()

    ##
    # @brief
    # @details
    # Reimplemented from MA.aligner.Module.execute.
    def execute(self, input_vec):
        sv_line = self.heap.pop()
        # @todo continue
        self.re_fill_heap()
        return sv_line
