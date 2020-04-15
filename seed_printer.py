from MS import *
from MA import *
from bokeh.plotting import figure, show
from bokeh_style_helper import *

class SeedPrinter(Module):
    def __init__(self, parameter_manager, name_a="Data", name_b="Ground Truth"):
        self.name_a = name_a
        self.name_b = name_b
        self.rendered = False

    ##
    # @brief Execute LineSweep for all given seeds.
    # @details
    # Reimplemented from MA.aligner.Module.execute.
    def execute(self, *input):
        assert(len(input) == 2)
        seeds_a = input[0]
        seeds_b = input[1]

        if not self.rendered:
            plot = figure(title="Seeds")
            def render(seed, name, dash=(10,0), width=2):
                if seed.on_forward_strand:
                    plot.line(x=[seed.start_ref, seed.start_ref + seed.size], y=[seed.start, seed.start + seed.size],
                              legend_label=name + " - forward", line_dash=dash, line_width=point_to_px(width))
                else:
                    plot.line(x=[seed.start_ref, seed.start_ref - seed.size], y=[seed.start, seed.start + seed.size],
                              legend_label=name + " - reverse", line_color="orange", line_dash=dash,
                              line_width=point_to_px(width))
            for seed in seeds_b:
                render(seed, self.name_b)
            for seed in seeds_a:
                render(seed, self.name_a, dash=(10,10), width=6)
            show(plot)

        self.rendered = True

        return None