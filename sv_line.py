from heap import *


class SvLine:
    def __init__(self, soc_id, ref_pos_start, ref_pos_end, query_pos, forw_strand, query_len):
        self.ref_pos_start = ref_pos_start
        self.ref_pos_end = ref_pos_end
        self.query_pos = query_pos
        self.nuc_seq_id, self.soc_index = soc_id
        self.forw_strand = forw_strand

        if not forw_strand and not query_pos is None:
            self.query_pos = query_len - self.query_pos

        self.num_supporting = None
        self.num_contradicting = None
        self.supp = None
        self.supporting = None

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
