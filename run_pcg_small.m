clear; clearvars; close all; beep off;
addpath('SuiteSparse-7.1.0/ssget');
addpath('utils');
warning('off', 'MATLAB:MKDIR:DirectoryExists');
rng(4751);

% Globals
global matvec_count;

% bldp options
oversampling = 20;
config_breg.method = 'krylov_schur';
config_breg.estimate_largest_with_nystrom = 1;
config_breg.tol = 1e-10;
config_breg.maxit = 200;
config_breg.oversampling = oversampling;
subspace_slack = 80;

config_evd.method = 'evd';
config_nys.method = 'nystrom';
config_nys.oversampling = oversampling;
config_nys_indef.method = 'indefinite_nystrom';
config_nys_indef.oversampling = oversampling;

% PCG parameters
tol_pcg = 1e-08;
maxit_pcg = 120;

% Paths and files
base_path = 'RESULTS/small';
mkdir(base_path);

csv_path = fullfile(base_path, 'csv_files');
mkdir(csv_path);
csv_header = "n,label,ratio,r,res,iter,flag,ctime,stime,matvecs,ksflag,cond,div\n";
csv_format = "%d,%s,%d,%d,%.2e,%d,%d,%d,%d,%d,%d,%.2e,%.2e\n";
plotting = Plotting();

options_file = fopen(fullfile(base_path, "options.txt"), "w");
fprintf(options_file,'Bregman Krylov-Schur options:\n');
fprintf(options_file,'estimate_largest_with_nystrom = %d\n', config_breg.estimate_largest_with_nystrom);
fprintf(options_file,'tol = %.1e\n', config_breg.tol);
fprintf(options_file,'maxit = %d\n', config_breg.maxit);
fprintf(options_file,'oversampling = %d\n', oversampling);
fprintf(options_file,'subspace_slack = %d\n', subspace_slack);
fprintf(options_file,'tol_pcg = %.1e\n', tol_pcg);
fprintf(options_file,'maxit_pcg = %d\n', maxit_pcg);
fclose(options_file);

% ichol and preconditioner parameters
default_options.type = 'nofill';
default_options.droptol = 0;  % ignored if 'type' is 'nofill'
default_options.michol = 'off';
default_options.diagcomp = 0;
options(1) = default_options;

ntols = 0;
droptols = logspace(-8, 1, ntols);
for i=2:numel(droptols)
    options(i).type = 'ict';
    options(i).droptol = droptols(i);  % ignored if 'type' is 'nofill'
    options(i).michol = 'on';
    options(i).diagcomp = 0;
end
diagcomp = 0.01;

% SuiteSparse matrices
names = ["494_bus", "1138_bus", "bcsstk04", "bcsstk05", "mesh2e1"];
names = [names, "bcsstk08", "bcsstk22", "bcsstm07", "nos5", "lund_a"];
suitesparse_criteria.names = names;
ids = SuitesSparseHelper.get(suitesparse_criteria);

rank_percentages = [0.01 0.05 0.1];
ratio_step = 0.25;  % Split between approximating positive and negative eigs

% Begin simulations
time_start = tic;
for i = 1:length(ids)
    id = ids(i);
    Prob = ssget(id);  % Prob is a struct (matrix, name, meta-data, ...)
    S = Prob.A;        % A is a symmetric sparse matrix
    S_action = @ (x) S_action_fn(S, x);
    label = replace(Prob.name, "/", "_");
    path_matrix = fullfile(base_path, label);

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
        G = full(Q \ S / Q');
        G = (G + G') / 2;
        eigenvalues = flip(sort(real(eig(G))));
        IplusG = I + G;
        nearness_measures = @ (p) condition_number_and_divergence(p, IplusG); 

        tic
        [~, flag_nopc, ~, iter_nopc, resvec_nopc] = pcg(S, b, tol_pcg, maxit_pcg);
        stime_nopc = toc;
        relres_nopc = resvec_nopc / norm_b;

        tic
        [~, flag_ichol, ~, iter_ichol, resvec_ichol] = pcg(S, b, tol_pcg, maxit_pcg, Q, Q');
        stime_ichol = toc;
        relres_ichol = resvec_ichol / norm_b;

        % Compute nearness
        div_nopc = bldp.bregman_divergence(S, I);
        div_ichol = bldp.bregman_divergence(I + G, I);

        % Write CSV header, unpreconditioned and ichol runs
        csv_out = fopen(fullfile(csv_path, [label '_ichol=' num2str(j) '.csv']), 'w');
        fprintf(csv_out, csv_header);
        fprintf(csv_out, csv_format, n, "nopc", -1, -1, relres_nopc(end), iter_nopc, flag_nopc, -1, stime_nopc, 0, -1, condest(S), div_nopc);
        fprintf(csv_out, csv_format, n, "ichol", -1, -1, relres_ichol(end), iter_ichol, flag_ichol, ctime_ichol, stime_ichol, 0, -1, condest(IplusG), div_ichol);
        
        has_already_failed = zeros(1, 1 + length(rank_percentages));
        % Loop for ranks
        for ridx = flip(1:numel(rank_percentages))
            any_success = 0;
            r = max(floor(n * rank_percentages(ridx)), 2);
            sketching_matrix = randn(n, r + oversampling);
            config_breg.sketching_matrix = sketching_matrix;

            % Exact SVD preconditioner
            config_evd.r = r;
            p_svd = bldp.svd_preconditioner(Q, S, config_evd);
            tic
            [~, flag_svd, ~, iter_svd, resvec_svd] = pcg(S, b, tol_pcg, maxit_pcg, p_svd.action);
            stime_svd = toc;
            relres_svd = resvec_svd / norm_b;
            [cond_svd, div_svd] = nearness_measures(p_svd);
            fprintf(csv_out, csv_format, n, "evd_exact", -1, r, relres_svd(end), iter_svd, flag_svd, p_svd.ctime, stime_svd, 1, -1, cond_svd, div_svd);

            % Exact Bregman preconditioner
            p_breg_exact = bldp.bregman_preconditioner(Q, S, config_evd);
            [cond_breg_exact, div_breg_exact] = nearness_measures(p_breg_exact);
            tic
            [~, flag_breg_exact, ~, iter_breg_exact, resvec_breg_exact] = pcg(S, b, tol_pcg, maxit_pcg, p_breg_exact.action);
            stime_breg_exact = toc;
            relres_breg = resvec_breg_exact/norm_b;
            fprintf(csv_out, csv_format, n, "breg_exact", -1, r, relres_breg(end), iter_breg_exact, flag_breg_exact, p_breg_exact.ctime, stime_breg_exact, 1, -1, cond_breg_exact, div_breg_exact);
            
            % Nyström of large positive eigenvalues
            config_nys.sketching_matrix = sketching_matrix;
            p_nys = bldp.svd_preconditioner(Q, S, config_nys);
            [cond_nys, div_nys] = nearness_measures(p_nys);
            tic
            [~, flag_nys, ~, iter_nys, resvec_nys] = pcg(S, b, tol_pcg, maxit_pcg, p_nys.action);
            stime_nys = toc;

            % Indefinite Nyström
            config_nys_indef.r = r;
            config_nys_indef.sketching_matrix = sketching_matrix;
            indefinite_nys_indef_fails = 0;
            if ~has_already_failed(1)
                p_nys_indef = bldp.svd_preconditioner(Q, S, config_nys_indef);
                try
                    [cond_nys_indef, div_nys_indef] = nearness_measures(p_nys_indef);
                catch
                    indefinite_nys_indef_fails = 1;
                end
            end
            if indefinite_nys_indef_fails || has_already_failed(1)
                p_nys_indef.c_time = -1;
                flag_nys_indef = 1; iter_nys_indef = -1; rv_nys_indef = -1;
                stime_nys_indef = -1; cond_nys_indef = -1; div_nys_indef = -1;
                has_already_failed(1) = 1;
            else
                tic
                [~, flag_nys_indef, ~, iter_nys_indef, resvec_nys_indef] = pcg(S, b, tol_pcg, maxit_pcg, p_nys_indef.action);
                stime_nys_indef = toc;
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
                    try
                        [cond_breg, div_breg] = nearness_measures(p_breg);
                    catch
                        bregman_krylovschur_fails = 1;
                    end
                end
                if bregman_krylovschur_fails || has_already_failed(1+ridx)
                    p_breg.diagnostics.nc = -1;
                    p_breg.ctime = -1;
                    flag_breg = 1; iter_breg = -1; rv_breg = -1;
                    stime_breg = -1; cond_breg = -1; div_breg = -1;
                    has_already_failed(1+ridx) = 1;
                else
                    tic
                    [~, flag_breg, ~, iter_breg, resvec_breg] = pcg(S, b, tol_pcg, maxit_pcg, p_breg.action);
                    stime_breg = toc;
                    rv_breg = resvec_breg(end)/norm_b;
                end
                fprintf(csv_out, csv_format, n, "breg", ...
                    config_breg.ratio, p_breg.diagnostics.nc, rv_breg, ...
                    iter_breg, flag_breg, p_breg.ctime, stime_breg, ...
                    matvec_count, p_breg.diagnostics.ks_flag, ...
                    cond_breg, div_breg);
                any_success = any_success || ~flag_breg;
            end
            fprintf(csv_out, csv_format, n, "nys", -1, r, resvec_nys(end)/norm_b, iter_nys, flag_nys, p_nys.ctime, stime_nys, 1, -1, cond_nys, div_nys);
            fprintf(csv_out, csv_format, n, "nys_indef", -1, r, resvec_nys_indef(end)/norm_b, iter_nys_indef, flag_nys_indef, p_nys_indef.ctime, stime_nys_indef, 1, -1, cond_nys_indef, div_nys_indef);

            % Plot PCG results
            path = fullfile(path_matrix, ['n=', num2str(n), '_r=', num2str(r), '_ichol=', num2str(j)]);
            mkdir(path);
            figure('Visible', 'off');

            plot_resvec(relres_nopc, plotting.nopc);
            plot_resvec(relres_ichol, plotting.ichol);
            plot_resvec(relres_breg, plotting.breg);
            plot_resvec(relres_svd, plotting.svd);

            xlabel('Iteration number', 'Interpreter', 'latex', 'FontSize', plotting.font_size_axes);
            ylabel('Relative residual', 'Interpreter', 'latex', 'FontSize', plotting.font_size_axes);
            ax = gca;
            ax.FontSize = plotting.font_size_ticks;
            grid on;
            ldg = legend;
            set(ldg, 'Interpreter', 'latex', 'FontSize', plotting.font_size_legend);
            set(ldg, 'Location', 'northoutside', 'Orientation', 'horizontal');
            ldg.AutoUpdate = 'off';
            yline(tol_pcg, 'r--');
            exportgraphics(gcf, fullfile(path, 'pcg_convergence.pdf'));
            hold off;
            
            % Plot Bregman vs SVD curves
            for ylog = [0 1]
                curves_path = fullfile(path, ['semilogy=', num2str(ylog)]);
                bldp_plot.plot_bregman_curves(eigenvalues, diag(p_svd.D),  diag(p_breg_exact.D), curves_path, ylog);
                bldp_plot.plot_svd_curve(eigenvalues, diag(p_svd.D), diag(p_breg_exact.D), curves_path, ylog);
            end
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

function [condn, div] = condition_number_and_divergence(p, M)
    condn = condest(bldp.SMW(p.U, p.D, p.V', M));
    div = bldp.bregman_divergence(M, bldp.SMW(p.U, p.D, p.V', speye(size(p.U, 1))));
end