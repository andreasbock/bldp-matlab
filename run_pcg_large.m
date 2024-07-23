clear; clearvars; close all; beep off;
addpath('SuiteSparse/ssget');
addpath('utils');
warning('off', 'MATLAB:MKDIR:DirectoryExists');
rng(4751);

% Globals
global matvec_count;

% bldp options
oversampling = 20;
config_breg.method = 'krylov_schur';
config_breg.estimate_largest_with_nystrom = 0;
config_breg.tol = 1e-04;
config_breg.maxit = 50;
config_breg.oversampling = oversampling;
config_nys.method = 'nystrom';
config_nys.oversampling = oversampling;
config_nys_indef.method = 'indefinite_nystrom';
config_nys_indef.oversampling = oversampling;
subspace_slack = 120;

% for csv files
label_nopc = "$I$";
label_ichol = "\texttt{ichol}";
label_breg_apx = @(alpha) strcat("$\PrecondBregAlpha{", num2str(alpha), "}$");
label_nys = "$\PrecondNys$";
label_nys_indef = "$\PrecondNysIndef$";

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
csv_header = "label,r,ratio,res,iter,flag,ctime,stime,matvecs,ksflag,cond,div\n";
csv_format = "%s,%d,%s,%.2e,%d,%d,%.2e,%.2e,%d,%d,%.2e,%.2e\n";
plotting = Plotting();

% ichol and preconditioner parameters
default_diagcomp = 0;
retry_diagcomp = 1e+03;

default_options.type = 'nofill';
default_options.droptol = 0;  % ignored if 'type' is 'nofill'
default_options.diagcomp = default_diagcomp;
options(1) = default_options;

ntols = 1;
drop_tols = logspace(-4, -4, ntols);
for i=1:numel(drop_tols)
    options(i+1).type = 'ict';
    options(i+1).droptol = drop_tols(i);  % ignored if 'type' is 'nofill'
    options(i+1).diagcomp = 0;
end

% dump options
options_file = fopen(fullfile(base_path, "options.txt"), "w");
fprintf(options_file, 'estimate_largest_with_nystrom = %d\n', config_breg.estimate_largest_with_nystrom);
fprintf(options_file, 'tol = %.1e\n', config_breg.tol);
fprintf(options_file, 'maxit = %d\n', config_breg.maxit);
fprintf(options_file, 'oversampling = %d\n', oversampling);
fprintf(options_file, 'subspace_slack = %d\n', subspace_slack);
fprintf(options_file, 'tol_pcg = %.1e\n', tol_pcg);
fprintf(options_file, 'maxit_pcg = %d\n', maxit_pcg);
fprintf(options_file, 'drop_tols = %d\n', drop_tols);
fprintf(options_file, 'default_diagcomp = %d\n', default_diagcomp);
fprintf(options_file, 'retry_diagcomp = %d\n', retry_diagcomp);
fclose(options_file);

% SuiteSparse matrices
names = ["Pres_Poisson", "bmwcra_1", "cfd1", "smt", "apache1", "thermal1"];
names = [names, "crankseg_1", "crankseg_2", "bcsstk17", "bcsstk18", "consph"];

suitesparse_criteria.names = names;
ids = SuitesSparseHelper.get(suitesparse_criteria);

rank_percentages = [0.0025, 0.0075];
ratio_step = 0.25;  % Split between approximating positive and negative eigs

% Begin simulations
time_start = tic;
for i = 1:length(ids)
    id = ids(i);
    Prob = ssget(id);  % Prob is a struct (matrix, name, meta-data, ...)
    S = Prob.A;        % A is a symmetric sparse matrix
    S_action = @ (x) S_action_fn(S, x);
    label = replace(Prob.name, "/", "_");

    n = size(S, 1);
    b = randn(n, 1);
    norm_b = norm(b);
    I = speye(n);
    config_breg.v = randn(n, 1);

    % Compute incomplete Cholesky
    for j = 1:numel(options)
        opts_ichol = options(j);
        try
            tic
            Q = ichol(S, opts_ichol);
            ctime_ichol = toc;
        catch ME
            fprintf('\t ichol failed (diagcomp = %d): %s\n \t Retrying with diagcomp = %d\n', opts_ichol.diagcomp, ME.message, retry_diagcomp);
            opts_ichol.diagcomp = retry_diagcomp;
            try
                tic
                Q = ichol(S, opts_ichol);
                ctime_ichol = toc;
            catch ME2
                fprintf('\t\t ichol STILL failed, abandoning: %s\n', ME.message);
                continue
            end
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
                        num2str(opts_ichol.diagcomp)];
        csv_out = fopen(fullfile(csv_path, [label ichol_string '.csv']), 'w');
        fprintf(csv_out, csv_header);
        fprintf(csv_out, csv_format, label_nopc, -1, "-1", relres_nopc(end), iter_nopc, flag_nopc, -1, stime_nopc, 0, -1);
        fprintf(csv_out, csv_format, label_ichol, -1, "-1", relres_ichol(end), iter_ichol, flag_ichol, ctime_ichol, stime_ichol, 0, -1);

        has_already_failed = zeros(1, 1 + length(rank_percentages));
        % Loop for ranks
        for ridx = flip(1:numel(rank_percentages))
            any_success = 0;
            r = max(floor(n * rank_percentages(ridx)), 2);
            sketching_matrix = randn(n, r + oversampling);
            config_breg.sketching_matrix = sketching_matrix;

            % Nyström of large positive eigenvalues
            config_nys.sketching_matrix = sketching_matrix;
            p_nys = bldp.svd_preconditioner(Q, S, config_nys);
            tic
            [~, flag_nys, ~, iter_nys, resvec_nys] = pcg(S, b, tol_pcg, maxit_pcg, p_nys.action);
            stime_nys = toc;

            % Indefinite Nyström
            config_nys_indef.r = r;
            config_nys_indef.sketching_matrix = sketching_matrix;
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
            any_success = any_success || ~flag_nys_indef ||  ~flag_nys;
            % Bregman
            for ratio = 0:ratio_step:1
                fprintf('[id = %s] %s, n = %d, r = %d, ratio = %f\n', num2str(id), label, n, r, ratio);
                matvec_count = 0;
                config_breg.ratio = ratio;
                config_breg.r = r;
                config_breg.subspace_dimension = r + subspace_slack;

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
                    p_breg.diagnostics.nc, num2str(ratio), rv_breg, ...
                    iter_breg, flag_breg, p_breg.ctime, stime_breg, ...
                    matvec_count, p_breg.diagnostics.ks_flag);
                any_success = any_success || ~flag_breg;
            end
            fprintf(csv_out, csv_format, label_nys, r, "-1", resvec_nys(end)/norm_b, iter_nys, flag_nys, p_nys.ctime, stime_nys, 1, -1);
            fprintf(csv_out, csv_format, label_nys_indef, r, "-1", resvec_nys_indef(end)/norm_b, iter_nys_indef, flag_nys_indef, p_nys_indef.ctime, stime_nys_indef, 1, -1);
        end
        fclose(csv_out);
    end
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