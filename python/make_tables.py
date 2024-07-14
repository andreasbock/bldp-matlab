import numpy as np
import pandas as pd

from pathlib import Path


def load_csv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, skipinitialspace=True)
    return df


def process_csv(csv_path_in, csv_path_out, maxit, tol=1e-10):
    csv = load_csv(csv_path_in)
    csv = csv.sort_values(by=['n', 'r'])

    tsvd_better = csv[csv['itersvd'] < csv['iterbreg']]
    tsvd_better_and_doesnt_fail = tsvd_better[(csv.flagsvd == 0) & (csv.flagbreg == 0)]
    if len(tsvd_better_and_doesnt_fail) > 0:
        raise Exception("TSVD better")

    nrow, _ = csv.shape
    pcs = ['ichol', 'svd', 'breg', 'rbreg']
    pc_names = {
        'cond': ['nopc'] + pcs,
        'iter': ['nopc'] + pcs,
        'div': pcs,
    }
    formats = {
        'cond': '{:.1e}',
        'iter': '{:d}',
        'div': '{:.1e}',
    }
    for i in range(nrow):
        csv['switch'][i] = int(i % 4 == 0)
        for pc in pcs:
            itr = 'iter' + pc
            if [csv[itr][i]] == 0:
                csv[itr][i] = maxit

        for prefix in ['cond', 'div']:
            indices = [prefix + tp for tp in pc_names[prefix]]
            _min = csv.loc[i, indices].astype(float).min()
            for idx in indices:
                if np.isinf(csv.loc[i, idx]):
                    csv[idx][i] = "-"
                elif abs(csv.loc[i, idx] - _min) < 1e-14:
                    csv[idx][i] = "\\textbf{" + formats[prefix].format(csv[idx][i]) + "}"
                else:
                    csv[idx][i] = formats[prefix].format(csv[idx][i])

        for prefix in ['iter']:
            indices = [prefix + tp for tp in pc_names[prefix]]
            nonzero_iters = csv.loc[i, indices].where(csv.loc[i, indices] > 0)
            _min = nonzero_iters.astype(float).min()
            for idx in pc_names[prefix]:
                iteridx = 'iter' + idx
                if abs(csv.loc[i, iteridx]) < 1e-14 or abs(csv.loc[i, 'res' + idx]) > tol:
                    csv[iteridx][i] = "-"
                elif abs(csv.loc[i, iteridx] - _min) < 1e-14:
                    csv[iteridx][i] = "\\textbf{" + formats[prefix].format(csv[iteridx][i]) + "}"
                else:
                    csv[iteridx][i] = formats[prefix].format(csv[iteridx][i])

        #for pc in pcs:
        #    res = 'res' + pc
        #    if np.isreal(csv[res][i]) and (csv[res][i] > tol or csv[itr][i] == 0):
        #        csv['iter' + pc][i] = "-"

    csv.to_csv(csv_path_out, index=False, float_format='%.1e')


if __name__ == '__main__':
    tol = 1e-09
    maxit = 100
    csv_in = Path('../RESULTS/small/results.csv')
    csv_out = Path('../RESULTS/small/results_out.csv')
    process_csv(csv_in, csv_out, maxit=maxit, tol=tol)
