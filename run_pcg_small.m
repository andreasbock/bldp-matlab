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
config_breg.estimate_largest_with_nystrom = 0;
config_breg.tol = 1e-10;
config_breg.maxit = 200;
config_breg.oversampling = oversampling;
subspace_slack = 80;

config_evd.method = 'evd';
config_nys.method = 'nystrom';
config_nys.oversampling = oversampling;
config_nys_indef.method = 'indefinite_nystrom';
config_nys_indef.oversampling = oversampling;

% for csv files
label_nopc = "$I$";
label_ichol = "\texttt{ichol}";
label_evd = "$\PrecondSVD$";
label_breg = "$\PrecondBreg$";
label_rbreg = "$\PrecondReverseBreg$";
label_breg_apx = @(alpha) strcat("$\PrecondBregAlpha{", num2str(alpha), "}$");
label_nys = "$\PrecondNys$";
label_nys_indef = "$\PrecondNysIndef$";

% PCG parameters
tol_pcg = 1e-09;
maxit_pcg = 100;

% Paths and files
base_path = 'RESULTS/small';
csv_path_all = fullfile(base_path, 'results.csv');

mkdir(base_path);
base_path = fullfile(base_path, strcat('nystrom_largest_eigs=', num2str(config_breg.estimate_largest_with_nystrom)));
mkdir(base_path);

% csv file per problem instance
csv_path = fullfile(base_path, 'csv_files');
mkdir(csv_path);
csv_header = "label,r,ratio,res,iter,flag,ctime,stime,matvecs,ksflag,cond,div\n";
csv_format = "%s,%d,%s,%.2e,%d,%d,%.2e,%.2e,%d,%d,%.2e,%.2e\n";
% csv file for all
csv_out_all = fopen(csv_path_all,'w');
csv_header_all = "Name,n,r,resnopc,resichol,ressvd,resbreg,resrbreg,iternopc,iterichol,itersvd,iterbreg,iterrbreg,condnopc,condichol,condsvd,condbreg,condrbreg,divnopc,divichol,divsvd,divbreg,divrbreg,flagnopc,flagichol,flagsvd,flagbreg,flagrbreg\n";
csv_format_all = '%s,%d,%d,%.2e,%.2e,%.2e,%.2e,%.2e,%d,%d,%d,%d,%d,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%d,%d,%d,%d,%d\n';
fprintf(csv_out_all, csv_header_all);
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
retry_diagcomp = 1e+02;
default_diagcomp = 0;
default_options.type = 'nofill';
default_options.droptol = 0;  % ignored if 'type' is 'nofill'
default_options.michol = 'off';
default_options.diagcomp = default_diagcomp;
options(1) = default_options;

ntols = 0;
droptols = logspace(-8, 1, ntols);
for i=2:numel(droptols)
    options(i).type = 'ict';
    options(i).droptol = droptols(i);  % ignored if 'type' is 'nofill'
    options(i).michol = 'on';
    options(i).diagcomp = default_diagcomp;
end

% SuiteSparse matrices
names = ["494_bus", "1138_bus", "bcsstk04", "662_bus", "bcsstk05"];
names = [names, "bcsstk08", "bcsstk22", "bcsstm07", "nos5", "lund_a"];
names = [names, "illc1850", "mesh2e1", "p0201", "lp_bandm"];
names = [names, "lp_sctap1", "lp_sctap3", "l9"];

suitesparse_criteria.names = names;
ids = SuitesSparseHelper.get(suitesparse_criteria);

rank_percentages = [0.01 0.05 0.1 0.175];
ratio_step = 0.25;  % Split between approximating positive and negative eigs

% Begin simulations
time_start = tic;
for i = 1:length(ids)
    id = ids(i);
    Prob = ssget(id);  % Prob is a struct (matrix, name, meta-data, ...)
    S = Prob.A;
    % Check if this is a least-squares problem; symmetrise if so
    if size(Prob.A, 1) > size(Prob.A, 2)
        S = Prob.A' * Prob.A;
    elseif size(Prob.A, 1) < size(Prob.A, 2)
        S = Prob.A * Prob.A';  % A is a symmetric sparse matrix
    end
    cond_S = condest(S);
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
        G = Q \ S / Q' - eye(n);
        G = full((G + G') / 2);
        eigenvalues = flip(sort(real(eig(G))));
        IplusG = I + G;
        cond_ichol = condest(IplusG);
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
        div_ichol = bldp.bregman_divergence(IplusG, I);

        % Write CSV header, unpreconditioned and ichol runs
        ichol_string = ['ichol_type=', opts_ichol.type, '_droptol=', ...
                        num2str(opts_ichol.droptol), '_diagcomp=', ...
                        num2str(opts_ichol.diagcomp)];
        csv_out = fopen(fullfile(csv_path, [label '_n=' num2str(n) '_' ichol_string '.csv']), 'w');
        fprintf(csv_out, csv_header);
        fprintf(csv_out, csv_format, label_nopc, -1, "-1", relres_nopc(end), iter_nopc, flag_nopc, -1, stime_nopc, 0, -1, cond_S, div_nopc);
        fprintf(csv_out, csv_format, label_ichol, -1, "-1", relres_ichol(end), iter_ichol, flag_ichol, ctime_ichol, stime_ichol, 0, -1, cond_ichol, div_ichol);
        
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
            fprintf(csv_out, csv_format, label_evd, r, "-1", relres_svd(end), iter_svd, flag_svd, p_svd.ctime, stime_svd, 1, -1, cond_svd, div_svd);

            % Exact Bregman preconditioner
            p_breg_exact = bldp.bregman_preconditioner(Q, S, config_evd);
            [cond_breg_exact, div_breg_exact] = nearness_measures(p_breg_exact);
            tic
            [~, flag_breg_exact, ~, iter_breg_exact, resvec_breg_exact] = pcg(S, b, tol_pcg, maxit_pcg, p_breg_exact.action);
            stime_breg_exact = toc;
            relres_breg = resvec_breg_exact/norm_b;
            fprintf(csv_out, csv_format, label_breg, r, "-1", relres_breg(end), iter_breg_exact, flag_breg_exact, p_breg_exact.ctime, stime_breg_exact, 1, -1, cond_breg_exact, div_breg_exact);
            
            % Exact reverse Bregman preconditioner
            p_rbreg_exact = bldp.reverse_bregman_preconditioner(Q, S, config_evd);
            [cond_rbreg_exact, div_rbreg_exact] = nearness_measures(p_rbreg_exact);
            tic
            [~, flag_rbreg_exact, ~, iter_rbreg_exact, resvec_rbreg_exact] = pcg(S, b, tol_pcg, maxit_pcg, p_rbreg_exact.action);
            stime_rbreg_exact = toc;
            relres_rbreg = resvec_rbreg_exact/norm_b;
            fprintf(csv_out, csv_format, label_rbreg, r, "-1", relres_rbreg(end), iter_rbreg_exact, flag_rbreg_exact, p_rbreg_exact.ctime, stime_rbreg_exact, 1, -1, cond_rbreg_exact, div_rbreg_exact);
            
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
                    rv_nys_indef = resvec_nys_indef(end)/norm_b;
                catch 
                    indefinite_nys_indef_fails = 1;
                end
            end
            if indefinite_nys_indef_fails || has_already_failed(1)
                p_nys_indef.c_time = -1;
                flag_nys_indef = 1; iter_nys_indef = -1; rv_nys_indef = -1;
                stime_nys_indef = -1; cond_nys_indef = -1; div_nys_indef = -1;
                has_already_failed(1) = 1; rv_nys_indef = -1;
            else
                tic
                [~, flag_nys_indef, ~, iter_nys_indef, resvec_nys_indef] = pcg(S, b, tol_pcg, maxit_pcg, p_nys_indef.action);
                stime_nys_indef = toc;
            end
            any_success = any_success || ~flag_nys_indef ||  ~flag_nys;
            % Bregman
            for ratio = 0:ratio_step:1
                fprintf('[id = %s] %s, n = %d, r = %d, ratio = %f, %s\n', num2str(id), label, n, r, ratio, ichol_string);
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
                fprintf(csv_out, csv_format, label_breg_apx(ratio), ...
                    p_breg.diagnostics.nc, num2str(ratio), rv_breg, iter_breg, ...
                    flag_breg, p_breg.ctime, stime_breg, ...
                    matvec_count, p_breg.diagnostics.ks_flag, ...
                    cond_breg, div_breg);
                any_success = any_success || ~flag_breg;
            end
            fprintf(csv_out, csv_format, label_nys, r, "-1", resvec_nys(end)/norm_b, iter_nys, flag_nys, p_nys.ctime, stime_nys, 1, -1, cond_nys, div_nys);
            fprintf(csv_out, csv_format, label_nys_indef, r, "-1", rv_nys_indef, iter_nys_indef, flag_nys_indef, p_nys_indef.ctime, stime_nys_indef, 1, -1, cond_nys_indef, div_nys_indef);

            % Plot PCG results
            path = fullfile(path_matrix, ['n=', num2str(n), '_r=', num2str(r), '_', ichol_string]);
            mkdir(path);
            figure('Visible', 'off');

            plot_resvec(relres_nopc, plotting.nopc);
            plot_resvec(relres_ichol, plotting.ichol);
            plot_resvec(relres_breg, plotting.breg);
            plot_resvec(relres_rbreg, plotting.rbreg);
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
            G_r = p_svd.U * p_svd.D * p_svd.V';
            G_r = (G_r + G_r')/2;
            e_svd = eig(G_r);
            ir = abs(e_svd)>1e-11;
            e_svd = e_svd(ir);
            G_breg = p_breg_exact.U * p_breg_exact.D * p_breg_exact.V';
            G_breg = (G_breg + G_breg')/2;
            e_breg_exact = eig(G_breg);
            ir = abs(e_breg_exact)>1e-11;
            e_breg_exact = e_breg_exact(ir);
            for ylog = [0 1]
                curves_path = fullfile(path, ['semilogy=', num2str(ylog)]);
                bldp_plot.plot_bregman_curves(eigenvalues, e_svd,  e_breg_exact, curves_path, ylog);
                bldp_plot.plot_svd_curve(eigenvalues, e_svd, e_breg_exact, curves_path, ylog);
            end

            % Dump to CSV
            fprintf(csv_out_all, csv_format_all, ...
                Prob.name, ...
                n, ...
                r, ...
                relres_nopc(end), ...
                relres_ichol(end), ...
                relres_svd(end), ...
                relres_breg(end), ...
                relres_rbreg(end), ...
                iter_nopc, ...
                iter_ichol, ...
                iter_svd, ...
                iter_breg_exact, ...
                iter_rbreg_exact, ...
                cond_S, ...
                cond_ichol, ...
                cond_svd, ...
                cond_breg_exact, ...
                cond_rbreg_exact, ...
                div_nopc, ...
                div_ichol, ...
                div_svd, ...
                div_breg_exact, ...
                div_rbreg_exact, ...
                flag_nopc, ...
                flag_ichol, ...
                flag_svd, ...
                flag_breg_exact, ...
                flag_rbreg_exact ...
            );
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
    div = bldp.bregman_divergence(M, speye(size(p.U, 1)) + p.U*p.D*p.V');
end