import numpy as np
import pandas as pd

from pathlib import Path


def load_csv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, skipinitialspace=True)
    return df


def process_csv(csv_path_in, csv_path_out, maxit, nranks, tol=1e-10):
    csv = load_csv(csv_path_in)
    csv = csv.sort_values(by=['n', 'r'])

    tsvd_better = csv[csv['itersvd'] < csv['iterbreg']]
    tsvd_better_and_doesnt_fail = tsvd_better[(csv.flagsvd == 0) & (csv.flagbreg == 0)]
    if len(tsvd_better_and_doesnt_fail) > 0:
        raise Exception("TSVD better")

    nrow, _ = csv.shape
    pcs_lowrank = ['svd', 'breg', 'rbreg']
    pcs = ['ichol'] + pcs_lowrank
    pcs_nopc = ['nopc'] + pcs
    pc_names = {
        'cond': pcs_nopc,
        'iter': pcs_nopc,
        'div': pcs,
    }
    formats = {
        'cond': '{:.1e}',
        'iter': '{:d}',
        'div': '{:.1e}',
    }
    for i in range(nrow):
        csv['switch'][i] = int(i % nranks == 0)
        for pc in pcs_nopc:
            itr = 'iter' + pc
            if csv[itr][i] == 0 or csv['flag' + pc][i] != 0:
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
            csv['switch'][i] = int(i % nranks == 0)
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

        eq_header = {
            'BeqR': ['iterbreg', 'iterrbreg'],
            'BeqS': ['iterbreg', 'itersvd'],
            'ReqS': ['iterrbreg', 'itersvd'],
        }
        for eq in eq_header.keys():
            if csv[eq][i] == 1:
                for v in eq_header[eq]:
                    csv[v][i] = f"{csv[v][i]}$^\\dagger$"

    csv.to_csv(csv_path_out, index=False, float_format='%.1e')


def process_csv_individual_tex(csv_path_in, csv_path_out, tol):
    csv = load_csv(csv_path_in)
    csv = csv.drop(['flag', 'ksflag'], axis=1).replace(-1, "-")
    csv['ratio'] = csv['ratio'].apply(lambda x: '{:.2f}'.format(x) if isinstance(x, float) else x)

    # Defining a function to apply
    mask_converged = csv.res > tol
    csv.loc[mask_converged, 'iter'] = "-"

    csv_string = csv.to_csv(header=None, index=False, float_format='%.1e', lineterminator='\\\\')
    csv_string = csv_string.replace("$\CSVPrecondBregAlpha{0}$", "\multirow{5}{*}{$\CSVPrecondBregAlpha{\alpha}$}")
    for ratio in ["0.25", "0.5", "0.75", "1"]:
        csv_string = csv_string.replace(f"$\CSVPrecondBregAlpha{{{ratio}}}$", "")
    csv_string = csv_string.replace(",", " & ").replace("\"", "")
    with open(csv_path_out, "w") as f:
        f.write(csv_string)


if __name__ == '__main__':
    csv_in = Path('../RESULTS/small/results.csv')
    csv_out = Path('../RESULTS/small/results_out.csv')
    process_csv(csv_in, csv_out, maxit=100, nranks=3, tol=1e-09)

    hpc_root = 'RESULTS_HPC'
    hpc_root_processed = f"{hpc_root}_tex"
    Path(f"../{hpc_root_processed}").mkdir(exist_ok=True)
    csv_in_paths = Path(f'../{hpc_root}/large').glob('nys*/csv_files/*.csv')

    for csv_in_path in csv_in_paths:
        print(f"Processing {csv_in_path.name}...")
        csv_out_path = Path(str(csv_in_path).replace(hpc_root, hpc_root_processed))
        csv_out_path.parent.mkdir(parents=True, exist_ok=True)
        process_csv_individual_tex(csv_in_path, csv_out_path, tol=1e-07)
