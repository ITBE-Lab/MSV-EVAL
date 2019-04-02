class SvJump:
    def __init__(self, seed, read_id, q_len):
        self.q_pos = None
        self.forw_strand = seed.on_forward_strand
        self.q_pos = seed.start + seed.size - 1
        self.ref_pos = None # noop
        if self.forw_strand:
            self.ref_pos = seed.start_ref + seed.size - 1
        else:
            self.ref_pos = seed.start_ref - seed.size + 1
        self.read_id = read_id
        self.destinations = []

    def add_destination(self, seed):
        class Destination:
            def __init__(self, parent, ref_pos, q_distance, switch_strands):
                self.parent = parent # creates pointer cycle
                self.ref_pos = ref_pos
                self.q_distance = q_distance
                self.switch_strands = switch_strands
        r_to = seed.start_ref
        if r_to == self.ref_pos:
            return False
        if self.forw_strand:
            if seed.start < self.q_pos:
                return False
            self.destinations.append(Destination(self,
                r_to, seed.start - self.q_pos, seed.on_forward_strand != self.forw_strand))
        else:
            if seed.start > self.q_pos:
                return False
            self.destinations.append(Destination(self,
                r_to, self.q_pos - seed.start, seed.on_forward_strand != self.forw_strand))
        return True