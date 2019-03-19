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
        return min(filter(lambda x: not x.empty(), self.heaps.values()), key=lambda x: x.peek().ref_pos_end)

    def peek(self):
        return self.first().peek()

    def pop(self):
        return self.first().pop()

    def empty(self):
        for x in self.heaps.values():
            if not x.empty():
                return False
        return True

    def __len__(self):
        return sum([len(x) for x in self.heaps.values()])

    def num_contradicting(self, ref_pos_start):
        contradicting = set()
        for sv_line in self.heaps[SeedSvLine]:
            if sv_line.ref_pos_start <= ref_pos_start:
                contradicting.add(sv_line.nuc_seq_id)
        return len(contradicting)

    def iter_contradicting(self):
        for sv_line in self.heaps[SeedSvLine]:
            if sv_line.ref_pos_start <= self.peek().ref_pos_start:
                yield sv_line

    def iter_supporting(self):
        for sv_line in self.heaps[type(self.peek())]:
            yield sv_line


##
# @brief overlaps sv lines and filters them
#
class SvLineFilter(VolatileModule):
    def __init__(self, db_name, pack, fm_index, parameter_manager):
        VolatileModule.__init__(self)  # this line is crucial DON'T DELETE ME

        self.min_support = parameter_manager.get_selected().by_name("min sv support").get()
        self.min_support_ratio = parameter_manager.get_selected().by_name(
            "min sv support ratio").get()

        self.heap = OpenSvLineHeaps()

        self.helper = SeedsToSVLines(
            db_name, pack, fm_index, parameter_manager)

        #
        # define different filters -> maybe move this to individual classes?
        #
        def type_not_seed():
            return not isinstance(self.heap.peek(), SeedSvLine)

        def support_ratio_filter():
            return self.heap.peek().supp >= self.min_support_ratio

        def support_filter():
            return self.heap.peek().num_supporting >= self.min_support

        # create a list of filters
        self.sv_line_filters = [
            type_not_seed,
            support_ratio_filter,
            support_filter
        ]

        self.next_sv_line = None
        if not self.helper.is_finished():
            self.next_sv_line = self.helper.execute(None)
        self.to_next_accepted_line()

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

    def decorate_curr_sv_line(self):
        sv_line = self.heap.peek()
        sv_line.num_supporting = len(
            self.heap.heaps[type(sv_line)])
        sv_line.num_contradicting = self.heap.num_contradicting(
            sv_line.ref_pos_start)
        sv_line.supp = sv_line.num_supporting / \
            (sv_line.num_contradicting + sv_line.num_supporting)
        sv_line.supporting = []
        for supporting_line in self.heap.iter_supporting():
            sv_line.supporting.append( (supporting_line.nuc_seq_id, supporting_line.query_pos) )

    # extract SV lines until the first one on the heap makes it through the filtering
    def to_next_accepted_line(self):
        self.re_fill_heap()
        while not self.is_finished():
            self.decorate_curr_sv_line()
            # print("to_next_accepted_line::while")
            filter_failed = False  # no filter has failed yet...
            for i, sv_line_filter in enumerate(self.sv_line_filters):
                if not sv_line_filter():
                    #print("Filter", i, "failed")
                    filter_failed = True  # continue to extract in the outer while loop
                    self.heap.pop()
                    break
            if filter_failed:
                #print("Filter failed")
                self.re_fill_heap()
                continue
            else:
                #print("all Filters passed")
                # if this point is reached and no filter returned False we exit the loop
                break

    ##
    # @brief
    # @details
    # Reimplemented from MA.aligner.Module.execute.
    def execute(self, input_vec):
        sv_line = self.heap.pop()
        self.to_next_accepted_line()
        return sv_line
        
##
# @brief overlaps sv lines and filters them
#
class SvLineOverlapper(VolatileModule):
    def __init__(self, db_name, pack, fm_index, parameter_manager):
        VolatileModule.__init__(self)  # this line is crucial DON'T DELETE ME

        self.helper = SvLineFilter(
            db_name, pack, fm_index, parameter_manager)

        self.last_sv_line = {
            SeedSvLine: None,
            GapStartSvLine: None,
            GapEndSvLine: None,
            OverlapSvLine: None
        }
        if self.helper.is_finished():
            self.set_finished()

    ##
    # @brief
    # @details
    # Reimplemented from MA.aligner.Module.execute.

    def execute(self, input_vec):
        while True:
            if self.helper.is_finished():
                remaining_sv_lines = list(filter(lambda x: not x is None, self.last_sv_line.values()))
                new_sv_line = remaining_sv_lines[0]
                self.last_sv_line[type(new_sv_line)] = None
                if len(remaining_sv_lines) == 1: # set_finished when returning the last sv line
                    self.set_finished()
                return new_sv_line
            else:
                new_sv_line = self.helper.execute(None)

                sv_line = self.last_sv_line[type(new_sv_line)]

                # if this is the first sv line of this type
                if sv_line is None:
                    self.last_sv_line[type(new_sv_line)] = new_sv_line
                    continue
                # check for non overlapping SV lines
                if sv_line.ref_pos_end < new_sv_line.ref_pos_start:
                    self.last_sv_line[type(new_sv_line)] = new_sv_line
                    break

                # overwrite overlapping sv line if necessary
                if new_sv_line.supp > sv_line.supp:
                    self.last_sv_line[type(new_sv_line)] = new_sv_line
                    continue
                if new_sv_line.supp == sv_line.supp and new_sv_line.num_supporting > sv_line.num_supporting:
                    self.last_sv_line[type(new_sv_line)] = new_sv_line
                    continue
                # extend the end position of the sv line otherwise
                sv_line.ref_pos_end = new_sv_line.ref_pos_end
                # write back the old sv line
                self.last_sv_line[type(new_sv_line)] = sv_line
        return sv_line


def add_sv_overlap_params(parameter_manager):
    parameter_manager.get_selected().register_parameter(libMA.AlignerParameterInt(
        "min sv support", "desc", 6, "Structural Variants Caller", 5))
    parameter_manager.get_selected().register_parameter(libMA.AlignerParameterDouble(
        "min sv support ratio", "desc", 6, "Structural Variants Caller", 1))


def test_SvLineOverlapper(fm_index, pack):
    print("testing SvLineOverlapper...")
    parameter_manager = ParameterSetManager()
    add_sv_line_params(parameter_manager)
    add_sv_overlap_params(parameter_manager)

    x = SvLineOverlapper("/MAdata/databases/sv_simulated",
                         pack, fm_index, parameter_manager)

    heap_size = 0
    num_steps = 0
    while not x.is_finished():
        sv_line = x.execute(None)
        print(sv_line.ref_pos_start, sv_line.ref_pos_end,
              sv_line.supp, sv_line.num_supporting, type(sv_line))
        num_steps += 1
        heap_size += len(x.helper.heap)
    print("success! Average heap size:", heap_size/num_steps)


if __name__ == "__main__":
    print("loading index...")
    fm_index = FMIndex()
    fm_index.load("/MAdata/genome/human/GRCh38.p12/ma/genome")
    pack = Pack()
    pack.load("/MAdata/genome/human/GRCh38.p12/ma/genome")
    print("done")

    test_SvLineOverlapper(fm_index, pack)

    print("done")
