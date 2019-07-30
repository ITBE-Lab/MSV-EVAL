from printColumns import print_columns
import math

class AnalyzeRuntimes:
    def __init__(self):
        self.times = {}

    def register(self, name, pledge, func=lambda x: x.exec_time):
        if not name in self.times:
            self.times[name] = []
        self.times[name].append( (pledge, func) )

    def analyze(self):
        print("runtime analysis:")
        data = []
        max_before_dot = 0
        for name, pledges in self.times.items():
            seconds = round(sum(func(pledge) for pledge, func in pledges), 3)
            max_before_dot = max(max_before_dot, int(math.log10(seconds)) )
            data.append([name, seconds])
        data = [(x, str(" "*int(max_before_dot-int(math.log10(y)))) + str(y) ) for x,y in data]
        data.insert(0, ["Module name", "runtime [s]"])
        print_columns(data)