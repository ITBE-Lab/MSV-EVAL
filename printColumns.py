def print_columns(data):
    col_width = [max([len(data[j][i]) for j in range(len(data))]) for i in range(len(data[0]))]
    first = True
    last_row = col_width
    for row in data:
        #if not first and cat != row[0]:
        #    print("| " + "".join(" "*l + " | " for l in col_width))
        #    cat = row[0]
        print("| " + "".join(
                (word.ljust(col_width[i]) if last_row[:i+1] != row[:i+1] \
                                          else " "*col_width[i]) + " | " for i, word in enumerate(row)
              ) )
        if first:
            print("-" * (sum(col_width) + len(col_width)*3 + 1))
            first = False
        last_row = row