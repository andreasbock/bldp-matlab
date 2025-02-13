clear; clearvars; close all; beep off;
addpath('utils');
warning('off', 'MATLAB:MKDIR:DirectoryExists');
rng(4751);

% Globals
global matvec_count;

% bldp options
oversampling = 30;
config_breg.method = 'krylov_schur';
config_breg.estimate_largest_with_nystrom = str2num(getenv("NYSTROM"));
config_breg.tol = 1e-02;
config_breg.oversampling = oversampling;
config_nys.method = 'nystrom';
config_nys.oversampling = oversampling;
config_nys_indef.method = 'nystrom';
config_nys_indef.oversampling = 1.5;
config_svd.method = 'krylov_schur';
config_svd.tol = 1e-02;

% for csv files
label_nopc = "$I$";
label_svd = "\CSVSVDKS";
label_ichol = "\CSVICHOL";
label_breg_apx = @(alpha) strcat("$\CSVPrecondBregAlpha{", num2str(alpha), "}$");
label_nys = "$\CSVPrecondNys$";
label_nys_indef = "$\CSVPrecondNysIndef$";

% PCG parameters
tol_pcg = 1e-07;
maxit_pcg = 350;

% Paths and files
base_path = 'RESULTS/large';
mkdir(base_path);
base_path = fullfile(base_path, strcat('nystrom_largest_eigs=', num2str(config_breg.estimate_largest_with_nystrom)));
mkdir(base_path);

csv_path = fullfile(base_path, 'csv_files');
mkdir(csv_path);
csv_header = "label,r,ratio,res,iter,flag,ctime,stime,matvecs,ksflag\n";
csv_format = "%s,%d,%s,%.2e,%d,%d,%.2e,%.2e,%d,%d\n";
plotting = Plotting();

% ichol and preconditioner parameters
default_options.type = 'nofill';
default_options.droptol = 0;  % ignored if 'type' is 'nofill'
options(1) = default_options;

drop_tols = [1e-01];
%options(1).type = 'ict';
%options(1).droptol = drop_tols(1);  % ignored if 'type' is 'nofill'

diagcomps = [0];
diagcomps(end+1) = -1;  % set later dynamically!
ndiagcomp = numel(diagcomps);

% Krylov-Schur parameters
subspace_slacks = [60, 100];
subspace_slack_svds = [60, 100];
maxits = [60, 100];

% dump options
options_file = fopen(fullfile(base_path, "options.txt"), "w");
fprintf(options_file, 'estimate_largest_with_nystrom = %d\n', config_breg.estimate_largest_with_nystrom);
fprintf(options_file, 'config_svd.tol = %.1e\n', config_svd.tol);
fprintf(options_file, 'config_breg.tol = %.1e\n', config_breg.tol);
fprintf(options_file, 'oversampling = %d\n', oversampling);
fprintf(options_file, 'maxits = %d\n', maxits);
fprintf(options_file, 'subspace_slacks = %d\n', subspace_slacks);
fprintf(options_file, 'subspace_slack_svds = %d\n', subspace_slack_svds);
fprintf(options_file, 'tol_pcg = %.1e\n', tol_pcg);
fprintf(options_file, 'maxit_pcg = %d\n', maxit_pcg);
fprintf(options_file, 'drop_tols = %d\n', drop_tols);
fprintf(options_file, 'diagcomp = %d\n', diagcomps);
fclose(options_file);

% Fetch matrix
suitesparse_criteria.names = [string(getenv("MATRIX_NAME"))];
ids = SuitesSparseHelper.get(suitesparse_criteria);

rank_percentages = [0.0025, 0.0075];

ratio_step = 0.25;  % Split between approximating positive and negative eigs

% Begin simulations
Prob = ssget(ids(1));  % Prob is a struct (matrix, name, meta-data, ...)
S = Prob.A;        % A is a symmetric sparse matrix
n = min(size(S, 1), size(S, 2));
b = randn(n, 1);
% Check if this is a least-squares problem; symmetrise if so
if size(S, 1) > size(S, 2)
    S = S' * S;
    b = S'*b;
elseif size(S, 1) < size(S, 2)
    S = S * S';
    b = S * b;
end
S_action = @ (x) S_action_fn(S, x);
label = replace(Prob.name, "/", "_");

norm_b = norm(b);
I = speye(n);
v = randn(n, 1);
config_breg.v = v;
config_svd.v = v;

% Compute incomplete Cholesky
time_start = tic;
for j = 1:numel(options)
    opts_ichol = options(j);
    ichol_success = 0;
    for kk = 1:numel(diagcomps)
        if kk == ndiagcomp
            % don't want to compute this if we have to
            diagcomps(end) = max(sum(abs(S), 2)./diag(S))-2;
        end
        opts_ichol.diagcomp = diagcomps(kk);
        try
            tic
            Q = ichol(S, opts_ichol);
            ctime_ichol = toc;
            ichol_success = 1;
            break;
        catch ME
            fprintf('\t ichol failed with diagcomp = %d: %s\n', opts_ichol.diagcomp, ME.message);
        end
    end
    if ~ichol_success
        fprintf('\t\t ichol failed for all options, abandoning.\n');
        continue;
    end
    tic
    [~, flag_nopc, ~, iter_nopc, resvec_nopc] = pcg(S, b, tol_pcg, maxit_pcg);
    stime_nopc = toc;
    relres_nopc = resvec_nopc / norm_b;

    tic
    [~, flag_ichol, ~, iter_ichol, resvec_ichol] = pcg(S, b, tol_pcg, maxit_pcg, Q, Q');
    stime_ichol = toc;
    relres_ichol = resvec_ichol / norm_b;
    
    % Write CSV header, unpreconditioned and ichol runs
    ichol_string = ['_ichol_type=', opts_ichol.type, '_droptol=', ...
                    num2str(opts_ichol.droptol), '_diagcomp=', ...
                    num2str(opts_ichol.diagcomp), ...
                    '_nnzS=', num2str(nnz(S)), ...
                    '_nnzQ=', num2str(nnz(Q))];
    csv_out = fopen(fullfile(csv_path, [label ichol_string '.csv']), 'w');
    fprintf(csv_out, csv_header);
    fprintf(csv_out, csv_format, label_nopc, -1, "-1", relres_nopc(end), iter_nopc, flag_nopc, -1, stime_nopc, 0, -1);
    fprintf(csv_out, csv_format, label_ichol, -1, "-1", relres_ichol(end), iter_ichol, flag_ichol, ctime_ichol, stime_ichol, 0, -1);

    has_already_failed = zeros(1, 1 + length(rank_percentages));
    % Loop for ranks
    for ridx = 1:numel(rank_percentages)
        any_success = 0;
        r = max(floor(n * rank_percentages(ridx)), 2);
        cr = round(r*config_nys_indef.oversampling);
        r_max = max([r + config_breg.oversampling, r + config_nys.oversampling, cr]);
        sketching_matrix = randn(n, r_max);

        % Nyström of large positive eigenvalues
        config_nys.r = r;
        config_nys.sketching_matrix = sketching_matrix(:, 1:r + config_nys.oversampling);
        p_nys = bldp.svd_preconditioner(Q, S, config_nys);
        tic
        [~, flag_nys, ~, iter_nys, resvec_nys] = pcg(S, b, tol_pcg, maxit_pcg, p_nys.action);
        stime_nys = toc;

        % Indefinite Nyström
        config_nys_indef.r = r;
        config_nys_indef.sketching_matrix = sketching_matrix(:, 1:cr);
        indefinite_nys_indef_fails = 0;
        if ~has_already_failed(1)
            p_nys_indef = bldp.svd_preconditioner(Q, S, config_nys_indef);
            tic
            [~, flag_nys_indef, ~, iter_nys_indef, resvec_nys_indef] = pcg(S, b, tol_pcg, maxit_pcg, p_nys_indef.action);
            stime_nys_indef = toc;
        else
            p_nys_indef.c_time = -1; flag_nys_indef = 1; 
            iter_nys_indef = -1; rv_nys_indef = -1;
            stime_nys_indef = -1;
            has_already_failed(1) = 1;
        end

        % SVD-based preconditioner using Krylov-Schur
        matvec_count = 0;
        config_svd.r = r;
	config_svd.maxit = maxits(ridx);
        config_svd.subspace_dimension = r + subspace_slack_svds(ridx);
        p_svd = bldp.svd_preconditioner(Q, S_action, config_svd);
        tic
        [~, flag_svd, ~, iter_svd, resvec_svd] = pcg(S, b, tol_pcg, maxit_pcg, p_svd.action);
        stime_svd = toc;
        fprintf(csv_out, csv_format, label_svd, rank(p_svd.D), "-1", resvec_svd(end)/norm_b, iter_svd, flag_svd, p_svd.ctime, stime_svd, matvec_count, p_svd.diagnostics.ks_flag);

        any_success = any_success || ~flag_nys_indef ||  ~flag_nys;
        % Bregman
	config_breg.maxit = maxits(ridx);
	config_breg.r = r;
	config_breg.subspace_dimension = r + subspace_slacks(ridx);
        for ratio = 0:ratio_step:1
            fprintf('%s, n = %d, r = %d, ratio = %f\n', label, n, r, ratio);
            matvec_count = 0;
            config_breg.r_largest = floor(config_breg.r * ratio);
            config_breg.sketching_matrix = sketching_matrix(:, 1:config_breg.r_largest + config_breg.oversampling);

            bregman_krylovschur_fails = 0;
            if ~has_already_failed(1+ridx)
                p_breg = bldp.bregman_preconditioner(Q, S_action, config_breg);
                tic
                [~, flag_breg, ~, iter_breg, resvec_breg] = pcg(S, b, tol_pcg, maxit_pcg, p_breg.action);
                stime_breg = toc;
                rv_breg = resvec_breg(end)/norm_b;
            else
                p_breg.diagnostics.nc = -1;
                p_breg.ctime = -1;
                flag_breg = 1; iter_breg = -1; rv_breg = -1;
                stime_breg = -1;
                has_already_failed(1+ridx) = 1;
            end
            fprintf(csv_out, csv_format, label_breg_apx(ratio), ...
                rank(p_breg.D), num2str(ratio), rv_breg, ...
                iter_breg, flag_breg, p_breg.ctime, stime_breg, ...
                matvec_count, p_breg.diagnostics.ks_flag);
            any_success = any_success || ~flag_breg;
        end
        fprintf(csv_out, csv_format, label_nys, r, "-1", resvec_nys(end)/norm_b, iter_nys, flag_nys, p_nys.ctime, stime_nys, r, -1);
        fprintf(csv_out, csv_format, label_nys_indef, r, "-1", resvec_nys_indef(end)/norm_b, iter_nys_indef, flag_nys_indef, p_nys_indef.ctime, stime_nys_indef, cr, -1);
    end
    fclose(csv_out);
end
time_end = toc(time_start);

disp("All done. Files are located in:");
disp(base_path);
disp("Total time in minutes:");
disp(time_end / 60);

function y = S_action_fn(S, x)
    global matvec_count;
    matvec_count = matvec_count + 1;
    y = S * x;
end
