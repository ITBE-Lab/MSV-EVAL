from heap import *


class SvLine:
    def __init__(self, line_id, ref_pos_start, ref_pos_end, query_pos):
        self.ref_pos_start = ref_pos_start
        self.ref_pos_end = ref_pos_end
        self.query_pos = query_pos
        self.line_id = line_id


class SeedSvLine(SvLine):
    def __init__(self, *args):
        SvLine.__init__(self, *args)


class GapStartSvLine(SvLine):
    def __init__(self, *args):
        SvLine.__init__(self, *args)


class GapEndSvLine(SvLine):
    def __init__(self, *args):
        SvLine.__init__(self, *args)


class OverlapSvLine(SvLine):
    def __init__(self, *args):
        SvLine.__init__(self, *args)


class OpenSvLineHeaps:
    def __init__(self):
        def lt_comp(x, y): return x.ref_pos_end < y.ref_pos_end
        self.heaps = {
            SeedSvLine = Heap(lt_comp),
            GapStartSvLine = Heap(lt_comp),
            GapEndSvLine = Heap(lt_comp),
            OverlapSvLine = Heap(lt_comp)
        }

    def push(self, sv_line):
        self.heaps[type(sv_line)].push(sv_line)

    def first(self):
        return min(self.heaps.values(), key=lambda x: x.peek().ref_pos_end)

    def peek(self):
        return self.first().peek()

    def pop(self)
    return self.first().pop()

    def empty():
        for x in self.heaps.values():
            if not x.empty():
                return False
        return True
