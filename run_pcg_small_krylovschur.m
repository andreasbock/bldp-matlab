clear; clearvars; close all; beep off;
addpath('SuiteSparse-7.1.0/ssget');
addpath('utils');
warning('off', 'MATLAB:MKDIR:DirectoryExists');

% Globals
global matvec_count;

% bldp options
config_breg.method = 'krylov_schur';
config_breg.estimate_largest_with_nystrom = 1;
config_breg.tol = 1e-10;
config_breg.maxit = 100;
config_breg.oversampling = 6;
subspace_slack = 20;

config_svd.method = 'nystrom';
config_svd.oversampling = 20;

% PCG parameters
tol_pcg = 1e-10;
maxit_pcg = 350;

% Paths and files
base_path = 'RESULTS/small_krylovschur';
mkdir(base_path);
csv_header = "n,label,ratio,r,res,iter,flag,ctime,stime,matvecs,ksflag\n";
csv_format = "%d,%s,%d,%d,%.2e,%d,%d,%d,%d,%d,%d\n";
options_file = fopen(fullfile(base_path, "options.txt"), "w");
fprintf(options_file,'Bregman Krylov-Schur options:\n');
fprintf(options_file,'estimate_largest_with_nystrom = %d\n', config_breg.estimate_largest_with_nystrom);
fprintf(options_file,'tol = %.1e\n', config_breg.tol);
fprintf(options_file,'maxit = %d\n', config_breg.maxit);
fprintf(options_file,'oversampling = %d\n', config_breg.oversampling);
fprintf(options_file,'subspace_slack = %d\n', subspace_slack);
fprintf(options_file,'tol_pcg = %.1e\n', tol_pcg);
fprintf(options_file,'maxit_pcg = %d\n', maxit_pcg);
fclose(options_file);

% ichol and preconditioner parameters
retry_diagcomp = 100;
default_diagcomp = 1;
default_opts_ichol.type = 'nofill';
default_opts_ichol.droptol = 0;  % ignored if 'type' is 'nofill'
default_opts_ichol.michol = 'off';
default_opts_ichol.diagcomp = default_diagcomp;
options(1) = default_opts_ichol;

ntols = 1;
droptols = logspace(-5, -3, ntols);
for i=1:numel(droptols)
    options(i+1).type = 'ict';
    options(i+1).droptol = droptols(i);  % ignored if 'type' is 'nofill'
    options(i+1).michol = 'on';
    options(i+1).diagcomp = 0;
end

% SuiteSparse matrices
names = ["494_bus", "1138_bus", "bcsstk04", "bcsstk05", "mesh2e1"];
names = [names, "bcsstk08", "bcsstk22", "bcsstm07", "nos5", "lund_a"];
suitesparse_criteria.names = names;
ids = SuitesSparseHelper.get(suitesparse_criteria);

ranks = [0.0075, 0.01];
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
        tic
        [~, flag_ichol, ~, iter_ichol, resvec_ichol] = pcg(S, b, tol_pcg, maxit_pcg, Q, Q');
        stime_ichol = toc;

        % Write CSV header, unpreconditioned and ichol runs
        csv_out = fopen(fullfile(base_path, [label '_ichol=' num2str(j) '.csv']), 'w');
        fprintf(csv_out, csv_header);
        fprintf(csv_out, csv_format, n, "nopc", -1, -1, resvec_nopc(end)/norm_b, iter_nopc, flag_nopc, -1, stime_nopc, 0, -1);
        fprintf(csv_out, csv_format, n, "ichol", -1, -1, resvec_ichol(end)/norm_b, iter_ichol, flag_ichol, ctime_ichol, stime_ichol, 0, -1);
    
        has_already_failed = zeros(1, 1 + length(ranks));
        % Loop for ranks
        for ridx = flip(1:numel(ranks))
            any_success = 0;
            r = round(n * ranks(ridx));
            if r == 0
                continue;
            end
            Omega = randn(n, r + config_svd.oversampling);
            % Nyström
            config_svd.Omega = Omega;
            if ~has_already_failed(1)
                p_nys = bldp.svd_preconditioner(Q, S, r, config_svd);
                tic
                [~, flag_nys, ~, iter_nys, resvec_nys] = pcg(S, b, tol_pcg, maxit_pcg, p_nys.action);
                stime_nys = toc;
                rv_nys = resvec_nys(end)/norm_b;
            else
                p_nys.c_time = -1;
                flag_nys = 1; iter_nys = -1; rv_nys = -1; stime_nys = -1;
                has_already_failed(1) = 1;
            end
            any_success = any_success || ~flag_nys;
            % Bregman
            for ratio = 0:ratio_step:1-ratio_step
                fprintf('[id = %s] %s, n = %d, r = %d, ratio = %f\n', ...
                        num2str(id), label, n, r, ratio);

                matvec_count = 0;
                config_breg.Omega = Omega;
                config_breg.ratio = ratio;
                config_breg.v = randn(n, 1);
                config_breg.subspace_dim = r + subspace_slack;
                if ~has_already_failed(1+ridx)
                    p_breg = bldp.bregman_preconditioner(Q, S_action, r, config_breg);
                    tic
                    [~, flag_breg, ~, iter_breg, resvec_breg] = pcg(S, b, tol_pcg, maxit_pcg, p_breg.action);
                    stime_breg = toc;
                    rv_breg = resvec_breg(end)/norm_b;
                else
                    p_breg.diagnostics.nc = -1;
                    p_breg.ctime = -1;
                    flag_breg = 1; iter_breg = -1; rv_breg = -1; stime_breg = -1;
                    has_already_failed(1+ridx) = 1;
                end
                fprintf(csv_out, csv_format, n, "breg", ...
                    config_breg.ratio, p_breg.diagnostics.nc, rv_breg, ...
                    iter_breg, flag_breg, p_breg.ctime, stime_breg, ...
                    matvec_count, p_breg.diagnostics.ks_flag ...
                );
                any_success = any_success || ~flag_breg;
            end
            fprintf(csv_out, csv_format, n, "nys", -1, r, rv_nys, iter_nys, flag_nys, p_nys.ctime, stime_nys, 1, -1);
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