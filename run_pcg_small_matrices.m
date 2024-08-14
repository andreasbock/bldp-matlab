clear; clearvars; close all; beep off;
addpath('utils');
warning('off', 'MATLAB:MKDIR:DirectoryExists');
rng(4751);

% Globals
global matvec_count;

% bldp options
config_breg.method = 'krylov_schur';
config_breg.estimate_largest_with_nystrom = 0;
config_breg.tol = 1e-10;
config_breg.maxit = 200;
config_breg.oversampling = 10;
subspace_slack = 50;

config_evd.method = 'evd';
config_nys.method = 'nystrom';
config_nys.oversampling = 20;
config_nys_indef.method = 'indefinite_nystrom';
config_nys_indef.oversampling = 1.5;

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
tol_pcg = 1e-10;
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
csv_header_all = "Name,n,r,resnopc,resichol,ressvd,resbreg,resrbreg,iternopc,iterichol,itersvd,iterbreg,iterrbreg,condnopc,condichol,condsvd,condbreg,condrbreg,divnopc,divichol,divsvd,divbreg,divrbreg,flagnopc,flagichol,flagsvd,flagbreg,flagrbreg,switch,BeqR,BeqS,ReqS\n";
csv_format_all = '%s,%d,%d,%.2e,%.2e,%.2e,%.2e,%.2e,%d,%d,%d,%d,%d,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%d,%d,%d,%d,%d,%d,%d,%d,%d\n';
fprintf(csv_out_all, csv_header_all);
plotting = Plotting();
eq_tol = 1e-11;

options_file = fopen(fullfile(base_path, "options.txt"), "w");
fprintf(options_file,'Bregman Krylov-Schur options:\n');
fprintf(options_file,'estimate_largest_with_nystrom = %d\n', config_breg.estimate_largest_with_nystrom);
fprintf(options_file,'tol = %.1e\n', config_breg.tol);
fprintf(options_file,'maxit = %d\n', config_breg.maxit);
if config_breg.estimate_largest_with_nystrom
    fprintf(options_file,'Bregman approx oversampling = %d\n', config_breg.oversampling);
end
fprintf(options_file,'Nyström oversampling = %d\n', config_nys.oversampling);
fprintf(options_file,'Indefinite Nyström oversampling = %d\n', config_nys_indef.oversampling);
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
    options(i).michol = 'off';
    options(i).diagcomp = default_diagcomp;
end

% SuiteSparse matrices
names = ["494_bus", "1138_bus", "662_bus", "bcsstk05"];
names = [names, "bcsstk08", "bcsstk22", "lund_a"];
names = [names, "illc1850", "mesh2e1", "lp_bandm"];
names = [names, "lp_sctap1", "lp_sctap3", "l9"];

suitesparse_criteria.names = names;
ids = SuitesSparseHelper.get(suitesparse_criteria);

rank_percentages = [0.01 0.05 0.1];
ratio_step = 0.25;  % Split between approximating positive and negative eigs

% Begin simulations
time_start = tic;
for i = 1:length(ids)
    id = ids(i);
    Prob = ssget(id);  % Prob is a struct (matrix, name, meta-data, ...)
    S = Prob.A;
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
    cond_S = condest(S);
    S_action = @ (x) S_action_fn(S, x);

    label = replace(Prob.name, "/", "_");
    path_matrix = fullfile(base_path, label);

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
        has_already_failed = zeros(1, 2 + length(rank_percentages));

        % Loop for ranks
        for ridx = flip(1:numel(rank_percentages))
            any_success = 0;
            r = max(floor(n * rank_percentages(ridx)), 2);
            cr = round(r*config_nys_indef.oversampling);
            r_max = max([r + config_breg.oversampling, r + config_nys.oversampling, cr]);
            sketching_matrix = randn(n, r_max);
            config_breg.sketching_matrix = sketching_matrix(:, 1:r + config_breg.oversampling);

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
            fprintf(csv_out, csv_format, label_breg, r, "-1", relres_breg(end), iter_breg_exact, flag_breg_exact, p_breg_exact.ctime, stime_breg_exact, -1, -1, cond_breg_exact, div_breg_exact);
            
            % Exact reverse Bregman preconditioner
            p_rbreg_exact = bldp.reverse_bregman_preconditioner(Q, S, config_evd);
            [cond_rbreg_exact, div_rbreg_exact] = nearness_measures(p_rbreg_exact);
            tic
            [~, flag_rbreg_exact, ~, iter_rbreg_exact, resvec_rbreg_exact] = pcg(S, b, tol_pcg, maxit_pcg, p_rbreg_exact.action);
            stime_rbreg_exact = toc;
            relres_rbreg = resvec_rbreg_exact/norm_b;
            fprintf(csv_out, csv_format, label_rbreg, r, "-1", relres_rbreg(end), iter_rbreg_exact, flag_rbreg_exact, p_rbreg_exact.ctime, stime_rbreg_exact, -1, -1, cond_rbreg_exact, div_rbreg_exact);

            % Plain Nyström
            config_nys.sketching_matrix = sketching_matrix(:, 1:r + config_nys.oversampling);
            nys_fails = 0;
            if ~has_already_failed(1)
                p_nys = bldp.svd_preconditioner(Q, S, config_nys);
                try
                    [cond_nys, div_nys] = nearness_measures(p_nys);
                catch 
                    min(eig(p_nys.U*p_nys.D*p_nys.V'))
                    nys_fails = 1;
                end
            end
            if nys_fails || has_already_failed(1)
                p_nys_indef.c_time = -1;
                flag_nys = 1; iter_nys = -1; rv_nys = -1;
                stime_nys = -1; cond_nys = -1; div_nys = -1;
                has_already_failed(1) = 1;
            else
                tic
                [~, flag_nys, ~, iter_nys, resvec_nys] = pcg(S, b, tol_pcg, maxit_pcg, p_nys.action);
                stime_nys = toc;
                rv_nys = resvec_nys(end)/norm_b;
            end

            % Park-Nakatsukasa Nyström
            config_nys_indef.r = r;
            config_nys_indef.sketching_matrix = sketching_matrix(:, 1:cr);
            nys_indef_fails = 0;
            if ~has_already_failed(2)
                p_nys_indef = bldp.svd_preconditioner(Q, S, config_nys_indef);
                try
                    [cond_nys_indef, div_nys_indef] = nearness_measures(p_nys_indef);
                catch 
                    min(eig(p_nys_indef.U*p_nys_indef.D*p_nys_indef.V'))
                    nys_indef_fails = 1;
                end
            end
            if nys_indef_fails || has_already_failed(2)
                p_nys_indef.c_time = -1;
                flag_nys_indef = 1; iter_nys_indef = -1; rv_nys_indef = -1;
                stime_nys_indef = -1; cond_nys_indef = -1; div_nys_indef = -1;
                has_already_failed(2) = 1;
            else
                tic
                [~, flag_nys_indef, ~, iter_nys_indef, resvec_nys_indef] = pcg(S, b, tol_pcg, maxit_pcg, p_nys_indef.action);
                stime_nys_indef = toc;
                rv_nys_indef = resvec_nys_indef(end)/norm_b;
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
                if ~has_already_failed(2+ridx)
                    p_breg = bldp.bregman_preconditioner(Q, S_action, config_breg);
                    try
                        [cond_breg, div_breg] = nearness_measures(p_breg);
                    catch
                        bregman_krylovschur_fails = 1;
                    end
                end
                if bregman_krylovschur_fails || has_already_failed(2+ridx)
                    p_breg.diagnostics.nc = -1;
                    p_breg.ctime = -1;
                    flag_breg = 1; iter_breg = -1; rv_breg = -1;
                    stime_breg = -1; cond_breg = -1; div_breg = -1;
                    has_already_failed(2+ridx) = 1;
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
            fprintf(csv_out, csv_format, label_nys, r, "-1", rv_nys, iter_nys, flag_nys, p_nys.ctime, stime_nys, r + config_nys.oversampling, -1, cond_nys, div_nys);
            fprintf(csv_out, csv_format, label_nys_indef, r, "-1", rv_nys_indef, iter_nys_indef, flag_nys_indef, p_nys_indef.ctime, stime_nys_indef, r, -1, cond_nys_indef, div_nys_indef);

            path = fullfile(path_matrix, ['n=', num2str(n), '_r=', num2str(r), '_', ichol_string]);
            mkdir(path);

            % Plot Bregman vs SVD curves
            G_r = p_svd.U * p_svd.D * p_svd.V';
            G_r = (G_r + G_r')/2;
            e_svd = eig(G_r);
            ir = abs(e_svd)>1e-11;
            e_svd = e_svd(ir);
            
            G_rbreg = p_rbreg_exact.U * p_rbreg_exact.D * p_rbreg_exact.V';
            G_rbreg = (G_rbreg + G_rbreg')/2;
            e_rbreg_exact = eig(G_rbreg);
            ir = abs(e_rbreg_exact)>1e-11;
            e_rbreg_exact = e_rbreg_exact(ir);

            G_breg = p_breg_exact.U * p_breg_exact.D * p_breg_exact.V';
            G_breg = (G_breg + G_breg')/2;
            e_breg_exact = eig(G_breg);
            ir = abs(e_breg_exact)>1e-11;
            e_breg_exact = e_breg_exact(ir);
            for ylog = [0 1]
                curves_path = fullfile(path, ['semilogy=', num2str(ylog)]);
                bldp_plot.plot_bregman_curves(eigenvalues, e_svd, e_breg_exact, curves_path, ylog);
                bldp_plot.plot_svd_curve(eigenvalues, e_svd, e_breg_exact, curves_path, ylog);
            end

            breg_equal_rbreg = norm(G_breg - G_rbreg) < eq_tol;
            breg_equal_svd = norm(G_breg - G_r) < eq_tol;
            rbreg_equal_svd = norm(G_rbreg - G_r) < eq_tol;
            
            eigs_sorted = sort(eigenvalues, 'descend');
            gamma = @(x) x - log(1 + x);
            nu = @(x) 1./(1 + x) + log(1 + x) - 1;
            % Plot spectra
            line_width = 1.8;
            font_size = 16;
            x = 1:n;
            figure('Visible', 'off');
            plot(x, eigs_sorted, 'Color', 'black', 'DisplayName', '$\tilde E$', 'LineStyle', '-', 'LineWidth', line_width);
            xlabel('$n$', 'Interpreter', 'latex', 'FontSize', font_size);
            ylabel('Eigenvalue', 'Interpreter', 'latex', 'FontSize', font_size);
            ldg = legend;
            set(ldg, 'Interpreter', 'latex', 'FontSize', font_size);
            grid on;
            exportgraphics(gcf, fullfile(path, 'eigenvalues.pdf'));
            hold off;

            % Plot spectra under the image of nu/gamma
            eps = 1e-16;
            line_width = 1.1;
            font_size = 16;
            x = 1:n;
            figure('Visible', 'off');
            semilogy(x, max(eps, gamma(eigs_sorted)), 'Color', 'red', 'DisplayName', '$\gamma$', 'LineStyle', '-', 'LineWidth', line_width); hold on;
            semilogy(x, max(eps, nu(eigs_sorted)), 'Color', 'blue', 'DisplayName', '$\nu$', 'LineStyle', '-', 'LineWidth', line_width);
            xlabel('$n$', 'Interpreter', 'latex', 'FontSize', font_size);
            ylim([eps 1e+06])
            ldg = legend;
            set(ldg, 'Interpreter', 'latex', 'Location', 'northwest', 'FontSize', font_size);
            grid on;
            exportgraphics(gcf, fullfile(path, 'gamma_nu.pdf'));
            hold off;

            % Plot PCG results
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
                flag_rbreg_exact, ...
                0, ...
                breg_equal_rbreg, ...
                breg_equal_svd, ...
                rbreg_equal_svd ...
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