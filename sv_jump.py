import math

class SvJump:
    def __init__(self, seed, read_id, q_len, jump_from_start):
        self.forw_strand = seed.on_forward_strand
        self.jump_from_start = jump_from_start
        self.read_id = read_id
        self.destinations = []

        self.q_pos = None  # noop
        self.ref_pos = None  # noop
        if jump_from_start:
            self.q_pos = seed.start
            self.ref_pos = seed.start_ref
        else:
            if self.forw_strand:
                self.q_pos = seed.start + seed.size - 1
                self.ref_pos = seed.start_ref + seed.size - 1
            else:
                self.q_pos = seed.start + seed.size - 1
                # @note direction mirrored on reference
                self.ref_pos = seed.start_ref - seed.size + 1

    def add_destination(self, seed):
        class Destination:
            def __init__(self, parent, ref_pos, q_distance, switch_strands, case, fuzziness_from_dir, fuzziness_to_dir):
                self.read_id = parent.read_id
                self.x = parent.ref_pos
                self.y = ref_pos
                self.ref_from_start = parent.ref_pos
                self.ref_from_end = parent.ref_pos
                self.ref_to_start = ref_pos
                self.ref_to_end = ref_pos
                f_1 = 1
                f = min(int(math.pow(max(abs(parent.ref_pos - ref_pos), q_distance), 1.5)/1000) + 1, 1000)
                if fuzziness_from_dir == "left":
                    self.ref_from_start -= f
                    self.ref_from_end += f_1
                else:
                    self.ref_from_end += f
                    self.ref_from_start -= f_1
                if fuzziness_to_dir == "down":
                    self.ref_to_start -= f
                    self.ref_to_end += f_1
                else:
                    self.ref_to_end += f
                    self.ref_to_start -= f_1
                self.switch_strands = switch_strands
                self.fuzziness_from_dir = fuzziness_from_dir
                self.fuzziness_to_dir = fuzziness_to_dir
                self.case = case
                self.score = 0.08 / math.log(q_distance + 1.5)
                self.q_distance = q_distance

            def __str__(self):
                return str(self.x) + ", " + str(self.y) + ", " + str(self.score) + ", " + str(self.case)

            # required for sorting but actually irrelevant
            def __lt__(self, other):
                return False

        r_to = None
        q_to = None
        case = None
        #fuzziness_from_dir = "right" if self.forw_strand else "left"
        #fuzziness_to_dir = "down" if seed.on_forward_strand else "up"
        #if seed.on_forward_strand != self.forw_strand: # == switch strands
        #    fuzziness_to_dir = fuzziness_from_dir
        fuzziness_from_dir = "right"
        fuzziness_to_dir = "down"
        if seed.on_forward_strand != self.forw_strand:
            if self.forw_strand and self.jump_from_start:
                fuzziness_from_dir = "left"
            if not self.forw_strand and not self.jump_from_start:
                fuzziness_from_dir = "left"
            if seed.on_forward_strand and self.jump_from_start:
                fuzziness_to_dir = "up"
            if not seed.on_forward_strand and not self.jump_from_start:
                fuzziness_to_dir = "up"
        # cases (0,2) and (0,3) in ppt matrix
        if self.jump_from_start and not seed.on_forward_strand:
            q_to = seed.start + seed.size - 1
            r_to = seed.start_ref - seed.size + 1  # @note direction mirrored on reference
            case = "(0,2-3)"
        # case (1,2) in ppt matrix
        elif not self.forw_strand and self.jump_from_start and seed.on_forward_strand:
            q_to = seed.start + seed.size - 1
            r_to = seed.start_ref + seed.size - 1
            case = "(1,2)"
        # case (2,1) in ppt matrix
        elif self.forw_strand and not self.jump_from_start and not seed.on_forward_strand:
            q_to = seed.start
            r_to = seed.start_ref
            case = "(2,1)"
        # cases (3,0) and (3,1) in ppt matrix
        elif not self.jump_from_start and seed.on_forward_strand:
            q_to = seed.start
            r_to = seed.start_ref
            case = "(3,0-1)"
        # in all other cases: return
        else:
            return None

        # finally we need to distinguish wether we jump up or downwards
        # (so that we dont get a negative jump distance on the query)
        # row 0 & 1 in ppt matrix
        if not self.jump_from_start and q_to >= self.q_pos:
            return Destination(self, r_to, q_to - self.q_pos,
                               seed.on_forward_strand != self.forw_strand, case,
                               fuzziness_from_dir, fuzziness_to_dir)
        # row 2 & 3 in ppt matrix
        elif self.jump_from_start and q_to <= self.q_pos:
            return Destination(self, r_to, self.q_pos - q_to,
                               seed.on_forward_strand != self.forw_strand, case,
                               fuzziness_from_dir, fuzziness_to_dir)
        return None
