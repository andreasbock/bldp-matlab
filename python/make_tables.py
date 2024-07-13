import numpy as np
import pandas as pd

from pathlib import Path


def load_csv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, skipinitialspace=True)
    return df


def process_csv(csv_path_in, csv_path_out, tol=1e-10):
    csv = load_csv(csv_path_in)
    csv = csv.sort_values(by=['n', 'r'])

    tsvd_better = csv[csv['itersvd'] < csv['iterbreg']]
    if len(tsvd_better) > 0:
        raise Exception("TSVD better")

    nrow, _ = csv.shape
    pc_names = ['nopc', 'ichol', 'svd', 'breg', 'rbreg']
    for i in range(nrow):

        for prefix in ['cond', 'div']:
            indices = [prefix + tp for tp in pc_names]
            _min = csv.loc[i, indices].astype(float).min()
            for idx in indices:
                if np.isinf(csv.loc[i, idx]):
                    csv[idx][i] = "-"
                elif abs(csv.loc[i, idx] - _min) < 1e-14:
                    csv[idx][i] = "\\textbf{" + '{:.1e}'.format(csv[idx][i]) + "}"
                else:
                    csv[idx][i] = '{:.1e}'.format(csv[idx][i])

        for pc in pc_names:
            res = 'res' + pc
            itr = 'iter' + pc

            if [csv[itr][i]] == 0:
                csv[itr][i] = 100

            if csv[res][i] > tol:
                csv[itr][i] = "-"

    csv.to_csv(csv_path_out, index=False, float_format='%.1e')


if __name__ == '__main__':
    tol = 1e-09
    csv_in = Path('../RESULTS/small/results.csv')
    csv_out = Path('../RESULTS/small/results_out.csv')
    process_csv(csv_in, csv_out, tol=tol)
