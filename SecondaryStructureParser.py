import sys
import os
import collections


class SecondaryStructureParser:

    def __init__(self, filename, autostart=True):
        self.filename = filename
        if not os.path.isfile(filename):
            self.filename = 'reading string directly'
            lines = filename
        else:
            self.filename = filename
            lines = None
        self.file_format = None
        self.parsed = {}
        self.valid_predictions = None
        if autostart:
            self.parse(lines=lines)
            self.validate()

    def guess_format(self, txt):
        if txt.startswith('#DeepConCNF_SS3'):
            return 'DeepConCNF_SS3'
        if txt.startswith('#DeepConCNF_SS8'):
            return 'DeepConCNF_SS8'
        raise ValueError('cannot guess type for file {}'.format(self.filename))

    def _set_valid_prediction(self):
        if self.file_format == 'DeepConCNF_SS3':
            self.valid_predictions = ('C', 'H', 'E')
        # TODO
        elif self.file_format == 'DeepConCNF_SS8':
            self.valid_predictions = ('C', 'H', 'E', 'L', 'T', 'S', 'G', 'B')

    def parse(self, lines=None):
        if lines is None:
            with open(self.filename, 'r') as f:
                lines = f.readlines()
        if len(lines) < 4:
            raise RuntimeError('File {} seems to be broken, it has less than 4 lines'.format(self.filename))
        self.file_format = self.guess_format(lines[0])
        self._set_valid_prediction()
        self._parser = getattr(self, '_parser_{}'.format(self.file_format))

        self.parsed = collections.OrderedDict()
        if self.file_format.startswith('DeepConCNF'):
            for line in lines[3:]:
                if line.strip() == '':
                    continue
                resp = self._parser(line)
                assert resp['id'] not in self.parsed
                self.parsed[resp['id']] = {'aa': resp['aa'],
                                           'prediction': resp['prediction'],
                                           'probabilities': resp.get('probabilities', ())}
        return self.parsed

    def validate(self):

        # make sure the IDs are incrementing by 1
        ids = list(self.parsed.keys())
        for i, id_ in enumerate(ids[:-1]):
            assert id_ + 1 in self.parsed, 'did not find residue ID: {}'.format(id_)

        # make sure the probabilities sum to 1
        for k, v in self.parsed.items():
            s = sum(v['probabilities'])
            assert 0.997 <= s <= 1.004, 'invalid sum of probabilities {} for residue {}'.format(s, k)

        # make sure the predictions are within the valid values
        for k, v in self.parsed.items():
            pred = v['prediction']
            assert pred in self.valid_predictions, 'invalid prediction "{}" for residue {}'.format(pred, k)

    def calculate_statistics(self):

        rel_occurence = collections.defaultdict(int)
        for v in self.parsed.values():
            rel_occurence[v['prediction']] += 1

        total_len = len(self.parsed)
        for k in self.valid_predictions:
            rel_occurence[k] /= total_len

        return dict(rel_occurence)

    def _parser(self):
        pass

    def _parser_DeepConCNF_SS3(self, line):
        return self._parser_DeepConCNF(line, 3)

    def _parser_DeepConCNF_SS8(self, line):
        return self._parser_DeepConCNF(line, 8)

    def _parser_DeepConCNF(self, line, n_predictions):
        aa_id = int(line[0:4])
        aa = line[5]
        pred = line[7]
        prob = line[11:]
        assert len(prob) == 6 * n_predictions + 1
        probs = [float(prob[0 + n * 6:6 + n * 6]) for n in range(n_predictions)]
        return {'id': aa_id, 'aa': aa, 'prediction': pred, 'probabilities': probs}


if __name__ == '__main__':
    if len(sys.argv) > 1:
        p = SecondaryStructureParser(sys.argv[1])
        print(p.calculate_statistics())