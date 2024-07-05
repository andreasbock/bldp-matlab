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

diagcomp = 0.01;
rank_percentages = [0.02 0.1];

% PCG parameters
tol_pcg = 1e-05;
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
    b = randn(1, n)';

    % Loop for ranks
    for ridx = 1:numel(rank_percentages)
        r_percentage = rank_percentages(ridx);
        r = max(floor(n * r_percentage), 2);
        
        % Loop over ichol options
        for j = 1:numel(options)
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
            G = Q \ S / Q' - I;
            G = (G + G')/2;
            G = Q \ S / Q';
            G = (G + G')/2;
            [eV, E] = eig(full(G));
            eg = real(diag(E)) - 1;

            % Absolute value: SVD
            config_svd.method = 'evd';
            p_svd = bldp.svd_preconditioner(Q, S, r, config_svd);
            
            % Bregman %
            % full eigendecomposition
            config_breg.method = 'evd';
            p_breg_full = bldp.bregman_preconditioner(Q, S, r, config_breg);

            % Krylov-Schur
            config_breg.method = 'krylov_schur';
            config_breg.estimate_largest_with_nystrom = 0;
            config_breg.tol = 1e-10;
            config_breg.maxit = 50;
            config_breg.subspace_dim = r + 10;
            config_breg.v = randn(n, 1);
            config_breg.r_largest = sum(sign(diag(p_breg_full.D))==1);
            p_breg_krs = bldp.bregman_preconditioner(Q,  @(x) S*x, r, config_breg);

            % Preconditioned conjugate gradient method
            [~, ~, ~, it_svd] = pcg(S, b, tol_pcg, maxit_pcg, p_svd.action);
            [~, ~, ~, it_full] = pcg(S, b, tol_pcg, maxit_pcg, p_breg_full.action);
            [~, ~, ~, it_krs] = pcg(S, b, tol_pcg, maxit_pcg, p_breg_krs.action);

            if abs(it_full - it_krs) > 7 || it_full > it_svd
                error("Test failed.");
            else
                disp("Test passed.");
            end
        end
    end
end
disp("ALL TESTS PASS.");
toc
