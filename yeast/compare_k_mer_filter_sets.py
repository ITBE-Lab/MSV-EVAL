from bokeh.plotting import figure, show, reset_output, ColumnDataSource

def load_filter_set(name):
    s = {}
    with open(name, "r") as in_file:
        for line in in_file.readlines():
            seq, cnt = line.split("\t")
            s[seq] = int(cnt)
    return s

def compare(read_filter_set, gt_filter_set, start, end, step, gt_threshold=1):
    print("loading...")
    gt_filter = load_filter_set(gt_filter_set)
    read_filter = load_filter_set(read_filter_set)

    print("gt_filter_cnt...")
    gt_filter_cnt = 0
    for _, cnt in gt_filter.items():
        if cnt > gt_threshold:
            gt_filter_cnt += 1
    print("gt_filter_cnt =", gt_filter_cnt)

    x_pos = []
    cnt_match = []
    cnt_extra = []
    cnt_missing = []
    print("cnt_match...")
    for idx in range(start, end, step):
        print(idx)
        match = 0
        extra = 0
        for seq, cnt in read_filter.items():
            if cnt >= idx:
                if seq in gt_filter and gt_filter[seq] > gt_threshold:
                    match += 1
                else:
                    extra += 1
        x_pos.append(idx)
        cnt_match.append(match)
        cnt_extra.append(extra)
        cnt_missing.append(gt_filter_cnt - match)

    return x_pos, cnt_match, cnt_missing, cnt_extra

def plot(x_pos, cnt_match, cnt_missing, cnt_extra):
    print("rendering")
    plot = figure(title="filter match", plot_width=1000, plot_height=1000)
    plot.line(x=x_pos, y=cnt_match, legend_label="#match", color="green")
    plot.line(x=x_pos, y=cnt_missing, legend_label="#missing", color="red")
    plot.line(x=x_pos, y=cnt_extra, legend_label="#extra", color="blue")
    show(plot)

if __name__ == "__main__":
    plot(*compare("tmp_k_mers.txt", "tmp_k_mers_assembly.txt", 0, 801, 50))