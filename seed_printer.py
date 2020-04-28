from MS import *
from MA import *
from bokeh.plotting import figure, show
from bokeh_style_helper import *

class SeedPrinter(Module):
    def __init__(self, parameter_manager, name_a="Data", name_b="Ground Truth"):
        self.name_a = name_a
        self.name_b = name_b

    ##
    # @brief Execute LineSweep for all given seeds.
    # @details
    # Reimplemented from MA.aligner.Module.execute.
    def execute(self, *input):
        assert(len(input) >= 2)
        seeds_a = input[0]
        seeds_b = input[1]

        plot = figure(title="Seeds")
        if len(input) > 2:
            helper_ret = input[2]
            for rect in helper_ret.rectangles:
                plot.quad(left=rect.x_axis.start, right=rect.x_axis.start+rect.x_axis.size,
                          bottom=rect.y_axis.start, top=rect.y_axis.start+rect.y_axis.size, color="lightgrey",
                          alpha=0.2)

        def render(seed, name, dash=(10,0), width=2):
            print(name, seed.start, seed.start_ref, seed.size, "forw" if seed.on_forward_strand else "rev")
            if seed.on_forward_strand:
                plot.line(x=[seed.start_ref, seed.start_ref + seed.size], y=[seed.start, seed.start + seed.size],
                            legend_label=name + " - forward", line_dash=dash, line_width=point_to_px(width))
            else:
                plot.line(x=[seed.start_ref, seed.start_ref - seed.size], y=[seed.start, seed.start + seed.size],
                            legend_label=name + " - reverse", line_color="orange", line_dash=dash,
                            line_width=point_to_px(width))
        for seed in seeds_b:
            render(seed, self.name_b, dash=(10,10), width=6)
        for seed in seeds_a:
            render(seed, self.name_a)
        show(plot)


        return None