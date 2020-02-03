from bokeh.layouts import column
from bokeh.plotting import figure, show, reset_output, ColumnDataSource
from MA import *

def analyze_jump_distance_distribution(sv_db, run_ids=None, resolution=5, x=None, y=None, w=None, h=None,
                                       max_jumps=100000, hard_max=1000, hard_min=10):
    if run_ids is None:
        run_ids = sv_db.newest_unique_runs( 3 )

    plots = []
    for cnt, run_id in enumerate(run_ids):
        assert sv_db.run_exists(run_id)
        print("sweeping:", sv_db.get_run_name(run_id))
        sweeper = None
        if not None in [x, y, w, h]:
            sweeper = SortedSvJumpFromSql(sv_db, sv_db.get_run_jump_id(run_id), x, y, w, h)
        else:
            sweeper = SortedSvJumpFromSql(sv_db, sv_db.get_run_jump_id(run_id))

        data = {}
        cnt = 0
        while sweeper.has_next_start():
            jump = sweeper.get_next_start()
            dist = int( abs(jump.from_pos - jump.to_pos) / resolution) * resolution
            if dist > hard_max or dist < hard_min:
                continue
            if not dist in data:
                data[dist] = 1
            else:
                data[dist] += 1
            cnt += 1
            if cnt >= max_jumps:
                break
        
        xs = []
        ys = []
        for x,y in data.items():
            xs.append(x)
            ys.append(y/cnt)

        if len(xs) > 0:
            plot = figure(width=1700, title=sv_db.get_run_name(run_id))
            if not len(plots) == 0:
                plot.x_range = plots[0].x_range
            plot.vbar(x=xs, top=ys, width=resolution*3/4)
            plots.append(plot)
    show(column(plots))


if __name__ == "__main__":
    sv_db = SV_DB("small_test_1", "open")
    analyze_jump_distance_distribution(sv_db)