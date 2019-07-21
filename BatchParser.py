    import sys
import os
import collections
import numpy as np
from SecondaryStructureParser import SecondaryStructureParser


class BatchParser:

    def __init__(self, files, suffixes=('.ss3', '.ss8'), autostart=True):
        if isinstance(files, str):
            if os.path.isdir(files):
                files = [os.path.join(files, f) for f in os.listdir(files) if f.endswith(suffixes)]
        self.files = sorted(files)

        self.results = {}
        self.statistics = collections.OrderedDict()
        if autostart:
            self.parse()

    def parse(self):
        self.results = collections.OrderedDict()
        for file in self.files:
            p = SecondaryStructureParser(file, autostart=True)
            self.results[file] = p.parsed
            self.statistics[file] = p.calculate_statistics()

    def calculate_statistics(self):
        predictions = set()
        for res in self.results.values():
            for residue in res.values():
                predictions.add(residue['prediction'])

        self.statistics['summary'] = collections.defaultdict(list)
        for k, v in self.statistics.items():
            if k == 'summary':
                continue
            for kk, vv in v.items():
                self.statistics['summary'][kk].append(vv)

        for k, v in self.statistics['summary'].items():
            self.statistics['summary'][k] = {
                'average': sum(v) / len(v),
                'median': np.median(v),
                'max': max(v),
                'min': min(v),
                'stdev': np.std(v)
            }
        self.statistics['total number of results'] = len(self.results)
        self.statistics['predictions'] = collections.defaultdict(float)
        total_len = 0
        for res in self.results.values():
            total_len += len(res.values())
            for residue in res.values():
                self.statistics['predictions'][residue['prediction']] += 1

        for k in self.statistics['predictions']:
            self.statistics['predictions'][k] /= total_len

        return self.statistics


if __name__ == '__main__':
    if len(sys.argv) == 1:
        b = BatchParser(os.getcwd())
    else:
        if len(sys.argv) > 2:
            suffixes = sys.argv[2]
        b = BatchParser(sys.argv[1], suffixes=suffixes)
        s = b.calculate_statistics()
