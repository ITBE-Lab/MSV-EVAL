from heap import *


class SvLine:
    def __init__(self, soc_id, ref_pos_start, ref_pos_end, query_pos, forw_strand):
        self.ref_pos_start = ref_pos_start
        self.ref_pos_end = ref_pos_end
        self.query_pos = query_pos
        self.nuc_seq_id, self.soc_index = soc_id
        self.forw_strand = forw_strand

        self.num_supporting = None
        self.num_contradicting = None
        self.supp = None
        self.supporting = None

    def merge(self, other):
        assert not self.supporting is None
        assert not other.supporting is None

        if other.ref_pos_start < self.ref_pos_start:
            self.ref_pos_start = other.ref_pos_start
        if other.ref_pos_end > self.ref_pos_end:
            self.ref_pos_end = other.ref_pos_end
        self.num_supporting += other.num_supporting
        self.num_contradicting += other.num_contradicting
        self.supp = self.num_supporting / (self.num_supporting + self.num_contradicting)
        self.supporting.extend(other.supporting)

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
