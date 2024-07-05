clear; clearvars; close all; beep off;
addpath('SuiteSparse-7.1.0/ssget');
addpath('utils');
warning('off', 'MATLAB:MKDIR:DirectoryExists');

tol_test = 1e-05;  % tolerance for tests

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
rank_percentages = [0.05 0.1];

% PCG parameters
tol_pcg = 1e-10;
maxit_pcg = 100;

% SuiteSparse matrices
suitesparse_criteria.names = ["494_bus", "1138_bus", "bcsstk05", "bcsstk08"];
ids = SuitesSparseHelper.get(suitesparse_criteria);

% Begin simulations
tic
for i = 1:length(ids)
    id = ids(i);
    Prob = ssget(id);  % Prob is a struct (matrix, name, meta-data, ...)
    S = Prob.A;        % A is a symmetric sparse matrix
    n = size(S, 1);
    I = speye(n);

    % Loop for ranks
    for ridx = 1:numel(rank_percentages)
        r_percentage = rank_percentages(ridx);
        r = max(floor(n * r_percentage), 2);
        % Loop over ichol options
        for j = 1:numel(options)
            fprintf("%s - %d (ichol = %d)\n", Prob.name, r, j);
            o = options(j);
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

            % Compute G and SVD/Bregman truncations
            G = full(Q \ S / Q' - I);
            G = (G + G')/2;
            [V, E] = eig(G);
            [~, gidx] = sort(abs(diag(E)));
            gidx_r = gidx(end-r+1:end);
            G_r = V(:, gidx_r)*E(gidx_r,gidx_r)*V(:,gidx_r)';
            bregman_curve = @(x) x - log(1 + x);
            [~, gidx_bm] = sort(bregman_curve(diag(E)));
            gidx_bm_r = gidx_bm(end-r+1:end);
            V_bregman = V(:,gidx_bm_r);
            E_bregman = E(gidx_bm_r,gidx_bm_r);
            G_bregman = V_bregman * E_bregman * V_bregman';

            % Absolute value: SVD
            config_abs.method = 'evd';
            p_svd = bldp.svd_preconditioner(Q, S, r, config_abs);
            error_svd_exact = norm(G_r - p_svd.U * p_svd.D * p_svd.U');

            % Absolute value: Nyström (note that r == n!)
            config_abs.method = 'nystrom';
            config_abs.Omega = randn(n, n);
            p_nys = bldp.svd_preconditioner(Q, S, n, config_abs);
            G_nys = p_nys.U * p_nys.D * p_nys.V';
            error_svd_nys = norm(G - G_nys);

            % Diagnostics
            % SVD vs KS (r << n) and SVD vs Nyström (r == n)
            error_svd = error_svd_exact > tol_test || error_svd_nys > tol_test;

            % Bregman
            config_breg.method = 'evd';
            % full eigendecomposition
            p_breg = bldp.bregman_preconditioner(Q, S, r, config_breg);
            error_bregman_exact = norm(G_bregman - p_breg.U * p_breg.D * p_breg.U');

            % Krylov-Schur
            config_breg.method = 'krylov_schur';
            config_breg.v = randn(n, 1);
            config_breg.estimate_largest_with_nystrom = 0;
            config_breg.tol = 1e-10;
            config_breg.maxit = 150;
            config_breg.subspace_dim = 150;

            % Compute ratio from analytical solution
            config_breg.r_largest = sum(sign(diag(p_breg.D))==1);
            p_breg_krs = bldp.bregman_preconditioner(Q, @(x) S*x, r, config_breg);
            error_bregman_ks = norm(G_bregman - p_breg_krs.U * p_breg_krs.D * p_breg_krs.U');
            
            % Diagnostics
            % Exact Bregman vs EVD and KS Bregman vs EVD
            error_bregman = error_bregman_exact > tol_test || error_bregman_ks > tol_test;
            
            % Test eigenvalues and eigenvectors
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
