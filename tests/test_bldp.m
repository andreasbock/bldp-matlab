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
tol_test = 1e-05;  % tolerance for tests

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
rank_percentages = [0.05 0.1];

%% PCG parameters
tol_pcg = 1e-10;
maxit_pcg = 100;

%% SuiteSparse matrices
names = ["494_bus", "1138_bus", "bcsstk05", "bcsstk08"];

suitesparse_criteria.n_max = 1200;
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

    %% Loop for ranks
    for ridx = 1:numel(rank_percentages)
        r_percentage = rank_percentages(ridx);
        r = max(floor(n * r_percentage), 2);
        
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

            G = Q \ S / Q';
            G = (G + G')/2;
            [eV, E] = eig(full(G));
            eg = real(diag(E)) - 1;

            %% Absolute value: SVD %%
            config_abs.evd = 1; config_abs.nystrom = 0;
            config_abs.krylov_schur = 0;
            p_svd = bldp.svd_preconditioner(Q, S, r, config_abs);
            %% Absolute value: Nyström
            config_abs.c = 1;
            config_abs.nystrom = 0; config_abs.nystrom = 1;
            config_abs.full_assembly = 1;
            % TODO!
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
            error_svd_ks = norm(p_ks_abs.V * p_ks_abs.D * p_ks_abs.V' - p_svd.V * p_svd.D * p_svd.V');
            % SVD vs Nyström (r == n)
            error_svd_nys = norm(p_nys.V * p_nys.D * p_nys.V' - G + I);

            error_svd = error_svd_ks > tol_test || error_svd_nys > tol_test;
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
            config_breg.bregman.ratio = 0.5;
            p_breg_krs = bldp.bregman_preconditioner(Q, S, r, config_breg);
            %% Diagnostics
            error_bregman_exact_eigs = norm(eg(p_breg_full.idx) - diag(p_breg_full.D));
            error_bregman_exact_evec = norm(eV(:,p_breg_full.idx) - p_breg_full.U);
            error_bregman_exact = max(error_bregman_exact_eigs, error_bregman_exact_evec);

            idx = flip([1:p_breg_krs.r_smallest, n-p_breg_krs.r_largest+1:n]);
            [p_breg_krs_eig, krs_idx] = sort(diag(p_breg_krs.D), 'descend');
            error_bregman_ks_eigs = norm(eg(idx) - p_breg_krs_eig);
            error_bregman_ks_evec = norm(abs(eV(:, idx)) - abs(p_breg_krs.U(:,krs_idx)));
            error_bregman_ks = max(error_bregman_exact_eigs, error_bregman_exact_evec);

            error_bregman = error_bregman_exact > tol_test || error_bregman_ks > tol_test;

            if error_svd || error_bregman
                error("Test failed.");
            else
                disp("Test passed.");
            end
        end
    end
end
disp("ALL TESTS PASS.");
toc
