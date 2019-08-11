import sys
import os
import collections
import pandas as pd


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
        if txt == '#\tAA\tSS\tHelix\tSheet\tCoil\n':
            return 'Porter3'
        if txt == '#\tAA\tSS\tG\tH\tI\tE\tB\tC\tS\tT\n':
            return 'Porter8'
        if '----- DISOPRED version 3' in txt:
            return 'Disopred3'
        raise ValueError('cannot guess type for file {}'.format(self.filename))

    def _set_valid_prediction(self):
        if self.file_format in ('DeepConCNF_SS3', 'Porter3'):
            self.valid_predictions = ('C', 'H', 'E')
        elif self.file_format in ('DeepConCNF_SS8', 'Porter8'):
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
            lines = lines[3:]
        elif self.file_format.startswith('Porter'):
            lines = lines[1:]
        elif self.file_format.startswith('Disopred'):
            lines = lines[3:]
        for line in lines:
            if line.strip() == '':
                continue
            resp = self._parser(line)
            assert resp['id'] not in self.parsed
            self.parsed[resp['id']] = {'aa': resp['aa']}
            if self.file_format.startswith('Disopred'):
                self.parsed[resp['id']] = {'disordered': resp['disordered'],
                                           'probability': resp['probability']}
            else:
                self.parsed[resp['id']] = {'prediction': resp['prediction'],
                                           'probabilities': resp.get('probabilities', ())}

        return self.parsed

    def validate(self):

        # make sure the IDs are incrementing by 1
        ids = list(self.parsed.keys())
        for i, id_ in enumerate(ids[:-1]):
            assert id_ + 1 in self.parsed, 'did not find residue ID: {}'.format(id_)

        # checks for ss3 and ss8 output
        # make sure the probabilities sum to 1
        if next(iter(self.parsed.values())).get('probabilities') is not None:
            for k, v in self.parsed.items():
                s = sum(v['probabilities'])
                assert 0.9969 <= s <= 1.004, 'invalid sum of probabilities {} for residue {}'.format(s, k)

            # make sure the predictions are within the valid values
            for k, v in self.parsed.items():
                pred = v['prediction']
                assert pred in self.valid_predictions, 'invalid prediction "{}" for residue {}'.format(pred, k)

        # checks for disopred output
        if next(iter(self.parsed.values())).get('disordered') is not None:
            for k, v in self.parsed.items():
                assert v['disordered'] in (True, False)
                assert 0 <= v['probability'] <= 1, v['probability']
                assert (v['probability'] > 0.5) == v['disordered'] or v['probability'] == 0.5, (self.filename, v['probability'], v['disordered'])

    def calculate_statistics(self):
        rel_occurence = collections.defaultdict(int)
        total_len = len(self.parsed)

        if self.file_format == 'Disopred3':
            for v in self.parsed.values():
                rel_occurence['disordered'] += int(v['disordered'])
            rel_occurence['disordered'] /= total_len
        else:
            for v in self.parsed.values():
                rel_occurence[v['prediction']] += 1

            for k in self.valid_predictions:
                rel_occurence[k] /= total_len

        return dict(rel_occurence)

    def to_df(self):
        return to_df(self.parsed)

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

    def _parser_Porter(self, line, n_predictions):
        cells = line.split('\t')
        assert len(cells) == 3 + n_predictions
        aa_id = int(cells[0])
        aa = cells[1]
        pred = cells[2]
        probs = [float(cell) for cell in cells[3:]]
        return {'id': aa_id, 'aa': aa, 'prediction': pred, 'probabilities': probs}

    def _parser_Porter3(self, line):
        return self._parser_Porter(line,3)

    def _parser_Porter8(self, line):
        return self._parser_Porter(line, 8)

    def _parser_Disopred3(self, line):
        aa_id = int(line[0:5])
        aa = line[6]
        disordered = line[8] == '*'
        probability = float(line[10:14])
        return {'id': aa_id, 'aa': aa, 'disordered': disordered, 'probability': probability}


def to_df(parsed):

    df = pd.DataFrame.from_dict(parsed, orient='index', columns=list(next(iter(parsed.values())).keys()))
    if 'probabilities' in df.columns:
        for i in range(len(df.head(1)['probabilities'].values[0])):
            df['probabilities_{}'.format(i)] = df['probabilities'].apply(lambda x: x[i])
        df.drop('probabilities', axis=1, inplace=True)

    return df


if __name__ == '__main__':
    if len(sys.argv) > 1:
        p = SecondaryStructureParser(sys.argv[1])
        print(p.calculate_statistics())

