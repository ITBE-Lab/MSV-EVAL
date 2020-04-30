from MS import *
from MA import *
from bokeh.plotting import figure, show
from ..bokeh_style_helper import *

class SeedPrinter(Module):
    def __init__(self, parameter_manager, name_a="Data", name_b="Ground Truth"):
        self.name_a = name_a
        self.name_b = name_b

    ##
    # @brief Execute LineSweep for all given seeds.
    # @details
    # Reimplemented from MA.aligner.Module.execute.
    def execute(self, *input):
        assert(len(input) >= 1)

        plot = figure(title="Seeds")
        if len(input) > 2:
            helper_ret = input[2]
            for rect in helper_ret.rectangles:
                plot.quad(left=rect.x_axis.start, right=rect.x_axis.start+rect.x_axis.size,
                          bottom=rect.y_axis.start, top=rect.y_axis.start+rect.y_axis.size, color="lightgrey",
                          alpha=0.2)

        def render(seeds, name, dash=(10,0), width=2):
            forw_x = []
            forw_y = []
            rev_x = []
            rev_y = []
            max_seed_size = max(seed.size for seed in seeds)
            for seed in seeds:
                if seed.size > max_seed_size*0.9:
                    print(name, seed.start, seed.start_ref, seed.size, "forw" if seed.on_forward_strand else "rev")
                if seed.on_forward_strand:
                    forw_x.append(seed.start_ref)
                    forw_x.append(seed.start_ref + seed.size)
                    forw_x.append(float("NaN"))
                    forw_y.append(seed.start)
                    forw_y.append(seed.start + seed.size)
                    forw_y.append(float("NaN"))
                else:
                    rev_x.append(seed.start_ref)
                    rev_x.append(seed.start_ref - seed.size)
                    rev_x.append(float("NaN"))
                    rev_y.append(seed.start)
                    rev_y.append(seed.start + seed.size)
                    rev_y.append(float("NaN"))
            plot.line(x=forw_x, y=forw_y,
                        legend_label=name + " - forward", line_dash=dash, line_width=point_to_px(width))
            plot.line(x=rev_x, y=rev_y,
                        legend_label=name + " - reverse", line_color="orange", line_dash=dash,
                        line_width=point_to_px(width))

        if len(input) > 1:
            seeds_b = input[1]
            render(seeds_b, self.name_b, dash=(10,10), width=6)

        seeds_a = input[0]
        render(seeds_a, self.name_a)

        show(plot)


        return None