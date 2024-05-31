clear; clearvars; close all; beep off;
addpath('SuiteSparse-7.1.0/ssget');
addpath('utils');
warning('off', 'MATLAB:MKDIR:DirectoryExists');

%% Paths and files
base_path = 'RESULTS';
mkdir(base_path);
csv_path = fullfile(base_path, 'results.csv');
csv_out = fopen(csv_path,'w');
fprintf(csv_out, "Name,n,r,resnopc,resichol,resscaled,resbreg,resreversebreg,iternopc,iterichol,iterscaled,iterbreg,iterreversebreg,condS,condichol,condscaled,condbreg,condreversebreg,divgichol,divgscaled,divgbreg,divgreversebreg,flagnopc,flagichol,flagscaled,flagbreg,flagreversebreg\n");
config = Config();

%% ichol and preconditioner parameters
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

%% PCG parameters
tol_pcg = 1e-10;
maxit_pcg = 100;

%% Krylov-Schur parameters
config.full_assembly = 1;
config.tol = 1e-05;
config.maxit = 40;
config.restart = 10;
%config.chunk_size = 4;

%% SuiteSparse matrices
names = ["494_bus"];  %, "1138_bus", "bcsstk04", "bcsstk05", "bcsstk18", "bcsstk08", "bcsstk22", "bcsstm07", "nos5", "lund_a", "mesh2e1", "bcsstk34"];
suitesparse_criteria.n_max = 5000;
suitesparse_criteria.n_min = 100;
suitesparse_criteria.symmetric = 1;
suitesparse_criteria.real = 1;
suitesparse_criteria.posdef = 1;
ids = suitesparse_helper.get(suitesparse_criteria, names);


%% Begin simulations
tic
for i = 1:length(ids)
    id = ids(i);
    Prob = ssget(id);  % Prob is a struct (matrix, name, meta-data, ...)
    S = Prob.A;        % A is a symmetric sparse matrix
    n = size(S, 1);
    I = speye(n);
    label = Prob.name;
    cond_S = condest(S);
    path_matrix = fullfile(base_path, label);

    %% Loop for ranks
    for ridx = 1:numel(rank_percentages)
        r_percentage = rank_percentages(ridx);
        r = max(floor(n * r_percentage), 2);
        subspace_iterations = floor((n / r));
        
        %% Loop over ichol options
        for j = 1:numel(options)
            o = options(j);
            opts_string = ['type=' num2str(o.type) '_droptol=' num2str(o.droptol) '_michol=' num2str(o.michol) '_diagcomp=', num2str(o.diagcomp)];
            fprintf('[id = %s] %s, n = %d, r = %d, (%s)\n', num2str(id), label, n, r, opts_string);
            path = fullfile(path_matrix, ['n=', num2str(n), '_r=', num2str(r), '_', opts_string]);

            %% Compute incomplete Cholesky
            try
                Q = ichol(S, o);
            catch ME
                fprintf('\t ichol failed: %s\n Retrying with diagcomp = %f\n', ME.message, diagcomp);
                o.diagcomp = diagcomp;
                path = fullfile(path_matrix, ['n=', num2str(n), '_r=', num2str(r), '_', opts_string]);
                try
                    Q = ichol(S, o);
                catch ME2
                    fprintf('\t\t ichol STILL failed, abandoning: %s\n', ME.message);
                end
                continue
            end
            G = Q \ S / Q' - I;
            [Gmin, Gmax] = bldp.extremal_eigenvalues(G, 5*r);
            if sign(Gmin) == sign(Gmax) || abs(Gmax + Gmin) < 1e-10
                fprintf('\t sign(Gmin) == sign(Gmax).\n');
                break
            end
            % compute ichol divergence and condition number
            if n < bldp.large_matrix_threshold()
                cond_ichol = condest(I + G);
                divg_ichol = bldp.bregman_divergence(I + G, I);
            else
                cond_ichol = 0;
                divg_ichol = 0;
            end

            %% Low-rank approximations %%            
            %% SVD (or Nyström)
            config_svd.nystrom = 1;
            config_svd.oversampling = 35;
            p_svd = bldp.svd_preconditioner(Q, S, r, config_svd);

            %% Bregman
            config.v = randn(n, 1);
            config.krylov_schur = 1;
            p_breg = bldp.bregman_preconditioner(Q, S, r, config);
            p_rbreg = bldp.reverse_bregman_preconditioner(Q, S, r, config);
            p_breg.diagnostics.ni

            %% Preconditioned conjugate gradient method
            b = randn(1, n)';
            norm_b = norm(b);
            pcg(S, b, tol_pcg, maxit_pcg, p_svd.action);
            pcg(S, b, tol_pcg, maxit_pcg, p_breg.action);

            [~, flag_nopc, ~, iter_nopc, resvec_nopc] = pcg(S, b, tol_pcg, maxit_pcg);
            [~, flag_ichol, ~, iter_ichol, resvec_ichol] = pcg(S, b, tol_pcg, maxit_pcg, Q, Q');
            [~, flag_svd, ~, iter_svd, resvec_svd] = pcg(S, b, tol_pcg, maxit_pcg, p_svd.action);
            [~, flag_breg, ~, iter_breg, resvec_breg] = pcg(S, b, tol_pcg, maxit_pcg, p_breg.action);
            [~, flag_rbreg, ~, iter_rbreg, resvec_rbreg] = pcg(S, b, tol_pcg, maxit_pcg, p_rbreg.action);

            relres_nopc = resvec_nopc / norm_b;
            relres_ichol = resvec_ichol / norm_b;
            relres_svd = resvec_svd / norm_b;
            relres_breg = resvec_breg / norm_b;
            relres_rbreg = resvec_rbreg / norm_b;
            mkdir(path);

            %% Plot PCG results
            figure('Visible', 'off');
            plot_resvec(relres_nopc, config.nopc);
            plot_resvec(relres_ichol, config.ichol);
            plot_resvec(relres_breg, config.breg);
            plot_resvec(relres_rbreg, config.rbreg);
            if n < bldp.large_matrix_threshold()
                config_absolute_value = config.svd;
            else
                config_absolute_value = config.nystrom;
            end
            plot_resvec(relres_svd, config_absolute_value);
            xlabel('Iteration number', 'Interpreter', 'latex', 'FontSize', config.font_size_axes);
            ylabel('Relative residual', 'Interpreter', 'latex', 'FontSize', config.font_size_axes);
            ax = gca;
            ax.FontSize = config.font_size_ticks;
            grid on;
            ldg = legend;
            set(ldg, 'Interpreter', 'latex', 'FontSize', config.font_size_legend);
            set(ldg, 'Location', 'northoutside', 'Orientation', 'horizontal');
            ldg.AutoUpdate = 'off';
            yline(tol_pcg, 'r--');
            exportgraphics(gcf, fullfile(path, 'pcg_convergence.pdf'));
            hold off;
            
            %% Condition numbers and divergence
            
            %% Plot Bregman vs SVD curves
            if n < bldp.large_matrix_threshold() && 0  % TODO
                eigenvalues = flip(sort(real(eig(full(G)))));
                for ylog = [0 1]
                    curves_path = fullfile(path, ['ylog = ', num2str(ylog)]);
                    bldp_plot.plot_bregman_curves(eigenvalues, p_svd.eig,  p_breg.eig, curves_path, ylog);
                    bldp_plot.plot_svd_curve(eigenvalues, p_svd.eig, p_breg.eig, curves_path, ylog);
                end
            end

            %% Dump to CSV
            fprintf(csv_out, '%s,%d,%d,%.2e,%.2e,%.2e,%.2e,%.2e,%d,%d,%d,%d,%d,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%.2e,%d,%d,%d,%d,%d\n', ...
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
                iter_breg, ...
                iter_rbreg, ...
                cond_S, ...
                cond_ichol, ...
                p_svd.nearness.cond, ...
                p_breg.nearness.cond, ...
                p_rbreg.nearness.cond, ...
                divg_ichol, ...
                p_svd.nearness.div, ...
                p_breg.nearness.div, ...
                p_rbreg.nearness.div, ...
                flag_nopc, ...
                flag_ichol, ...
                flag_svd, ...
                flag_breg, ...
                flag_rbreg ...
            );
        end
    end
end
fclose(csv_out);
toc
