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
config_breg = Config();
tol_test = 1e-09;  % tolerance for tests

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

%% SuiteSparse matrices
names = ["494_bus"];%, "1138_bus", "bcsstk04", "bcsstk05", "bcsstk18", "bcsstk08", "bcsstk22", "bcsstm07", "nos5", "lund_a", "mesh2e1", "bcsstk34"];

suitesparse_criteria.n_max = 1000;
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
            %% Compute incomplete Cholesky
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
            G = (G + G')/2;
            [Gmin, Gmax] = bldp.extremal_eigenvalues(G, 5*r);
            if sign(Gmin) == sign(Gmax) || abs(Gmax + Gmin) < 1e-10
                fprintf('\t sign(Gmin) == sign(Gmax).\n');
                break
            end
            %e = flip(sort(eig(G)));
            %% Absolute value: SVD %%
            config_abs.evd = 1; config_abs.nystrom = 0;
            config_abs.krylov_schur = 0;
            p_svd = bldp.svd_preconditioner(Q, S, r, config_abs);
            %% Absolute value: Nyström
            config_abs.c = 1;
            config_abs.nystrom = 0; config_abs.nystrom = 1;
            config_abs.full_assembly = 1;
            p_nys = bldp.svd_preconditioner(Q, S, n, config_abs);
            %% Absolute value: Krylov-Schur
            config_abs.tol = 1e-10;
            config_abs.maxit = 150;
            config_abs.restart = 150;
            config_abs.v = randn(n, 1);
            config_abs.krylov_schur = 1; config_abs.nystrom = 0; config_abs.nystrom = 0;
            p_ks_abs = bldp.svd_preconditioner(Q, S, r, config_abs);
            %% Diagnostics
            % SVD vs KS (r << n)
            error_svd_ks = norm(p_ks_abs.V * p_ks_abs.D * p_ks_abs.V' - p_svd.V * p_svd.D * p_svd.V')
            % SVD vs Nyström (r == n)
            error_svd_nys = norm(p_nys.V * p_nys.D * p_nys.V' - G)

            %% Bregman %%
            config_breg.tol = 1e-10;
            config_breg.maxit = 100;
            config_breg.restart = 100;
            config_breg.v = randn(n, 1);
            config_breg.krylov_schur = 0;
            % full eigendecomposition
            p_breg_full = bldp.bregman_preconditioner(Q, S, r, config_breg);
            % Krylov-Schur
            config_breg.krylov_schur = 1;
            config_breg.bregman.min_eps = 0.1;
            config_breg.bregman.max_eps = 0.1;
            config_breg.bregman.ratio = 0.8;
            p_breg_krs = bldp.bregman_preconditioner(Q, S, r, config_breg);
            %% Diagnostics
            p_breg_krs.diagnostics
            %error_bregman_ks = norm(p_breg_full.D - p_breg_krs.D)

            if error_svd_ks > tol_test || error_svd_nys > tol_test
                error("Test failed.");
            else
                fprintf("Test passed.\n");
            end
        end
    end
end
fclose(csv_out);
toc
