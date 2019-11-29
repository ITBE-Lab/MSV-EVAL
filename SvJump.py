class SvJump:
    def __init__(self, seed, read_id):
        self.ref_pos = seed.start_ref + seed.size
        self.q_pos = seed.start + seed.size
        self.forw_strand = seed.on_forward_strand
        self.read_id = read_id
        self.destinations = []

    def add_destination(self, seed):
        class Destination:
            def __init__(self, ref_pos, q_distance, switch_strands):
                self.ref_pos = ref_pos
                self.q_distance = q_distance
                self.switch_strands = switch_strands
        assert seed.start >= self.q_pos
        self.destinations.append(Destination(
            seed.start_ref, seed.start - self.q_pos, seed.on_forward_strand != self.forw_strand))