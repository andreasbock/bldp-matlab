clear; clearvars; close all; beep off;
addpath('SuiteSparse-7.1.0/ssget');
addpath('utils');
warning('off', 'MATLAB:MKDIR:DirectoryExists');

% Paths and files
base_path = 'RESULTS/small';
mkdir(base_path);
csv_path = fullfile(base_path, 'results.csv');
csv_out = fopen(csv_path,'w');
csv_header = "Name,n,r,resnopc,resichol,resscaled,resbreg,iternopc,iterichol,iterscaled,iterbreg,condS,condichol,condscaled,condbreg,divgichol,divgscaled,divgbreg,flagnopc,flagichol,flagscaled,flagbreg\n";
csv_format = '%s,%d,%d,%.2e,%.2e,%.2e,%.2e,%d,%d,%d,%d,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%d,%d,%d,%d\n';
fprintf(csv_out, csv_header);
plotting = Plotting();

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
rank_percentages = [0.01 0.05 0.1];

% PCG parameters
tol_pcg = 1e-7;
maxit_pcg = 50;

% SuiteSparse matrices
names = ["494_bus", "1138_bus", "bcsstk04", "bcsstk05", "mesh2e1"];
names = [names, "bcsstk08", "bcsstk22", "bcsstm07", "nos5", "lund_a"];
suitesparse_criteria.names = names;
ids = SuitesSparseHelper.get(suitesparse_criteria);

% bldp options
config_svd.method = 'evd';
config_breg.method = 'evd';

% Begin simulations
tic
for i = 1:length(ids)
    id = ids(i);
    Prob = ssget(id);  % Prob is a struct (matrix, name, meta-data, ...)
    S = Prob.A;        % A is a symmetric sparse matrix
    n = size(S, 1);
    I = speye(n);
    b = randn(1, n)';
    norm_b = norm(b);
    label = Prob.name;
    cond_S = condest(S);
    path_matrix = fullfile(base_path, label);

    % Loop over ichol options
    for j = 1:numel(options)
    
        o = options(j);
        %opts_string = ['type=' num2str(o.type) '_droptol=' num2str(o.droptol) '_michol=' num2str(o.michol) '_diagcomp=', num2str(o.diagcomp)];    
        % Compute incomplete Cholesky
        try
            Q = ichol(S, o);
        catch ME
            fprintf('\t ichol failed: %s\n Retrying with diagcomp = %f\n', ME.message, diagcomp);
            o.diagcomp = diagcomp;
            try
                Q = ichol(S, o);
            catch ME2
                fprintf('\t\t ichol STILL failed, abandoning: %s\n', ME.message);
            end
            continue
        end
        G = Q \ S / Q' - I;
        
        % compute ichol divergence and condition number
        cond_ichol = condest(I + G);
        divg_ichol = bldp.bregman_divergence(I + G, I);
        
        % Loop for ranks
        for ridx = 1:numel(rank_percentages)
            r_percentage = rank_percentages(ridx);
            r = max(floor(n * r_percentage), 2);
            fprintf('[id = %s] %s, n = %d, r = %d\n', num2str(id), label, n, r);

            subspace_iterations = floor((n / r));
            
            % Low-rank approximations
            % SVD
            p_svd = bldp.svd_preconditioner(Q, S, r, config_svd);
            % Bregman
            p_breg = bldp.bregman_preconditioner(Q, S, r, config_breg);

            % Preconditioned conjugate gradient method
            [~, flag_nopc, ~, iter_nopc, resvec_nopc] = pcg(S, b, tol_pcg, maxit_pcg);
            [~, flag_ichol, ~, iter_ichol, resvec_ichol] = pcg(S, b, tol_pcg, maxit_pcg, Q, Q');
            [~, flag_svd, ~, iter_svd, resvec_svd] = pcg(S, b, tol_pcg, maxit_pcg, p_svd.action);
            [~, flag_breg, ~, iter_breg, resvec_breg] = pcg(S, b, tol_pcg, maxit_pcg, p_breg.action);
            
            relres_nopc = resvec_nopc / norm_b;
            relres_ichol = resvec_ichol / norm_b;
            relres_svd = resvec_svd / norm_b;
            relres_breg = resvec_breg / norm_b;

            path = fullfile(path_matrix, ['n=', num2str(n), '_r=', num2str(r), '_ichol=', num2str(j)]);
            mkdir(path);


            % Plot PCG results
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
            G = full(G);
            G = (G + G') / 2;
            eigenvalues = flip(sort(real(eig(G))));
            for ylog = [0 1]
                curves_path = fullfile(path, ['semilogy=', num2str(ylog)]);
                bldp_plot.plot_bregman_curves(eigenvalues, diag(p_svd.D),  diag(p_breg.D), curves_path, ylog);
                bldp_plot.plot_svd_curve(eigenvalues, diag(p_svd.D), diag(p_breg.D), curves_path, ylog);
            end

            % Divergences
            divg_svd = bldp.bregman_divergence(I + G, bldp.SMW(p_svd.U, p_svd.D, p_svd.U', I));
            divg_breg = bldp.bregman_divergence(I + G, bldp.SMW(p_breg.U, p_breg.D, p_breg.U', I));
            
            % Dump to CSV
            fprintf(csv_out, csv_format, ...
                Prob.name, ...
                n, ...
                r, ...
                relres_nopc(end), ...
                relres_ichol(end), ...
                relres_svd(end), ...
                relres_breg(end), ...
                iter_nopc, ...
                iter_ichol, ...
                iter_svd, ...
                iter_breg, ...
                cond_S, ...
                cond_ichol, ...
                condest(bldp.SMW(p_svd.U, p_svd.D, p_svd.U', I + G)), ...
                condest(bldp.SMW(p_breg.U, p_svd.D, p_svd.U', I + G)), ...
                divg_ichol, ...
                divg_svd, ...
                divg_breg, ...
                flag_nopc, ...
                flag_ichol, ...
                flag_svd, ...
                flag_breg ...
            );
        end
    end
end
fclose(csv_out);
toc
