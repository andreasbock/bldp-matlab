clear; clearvars; close all; beep off;

% for suitesparse_helper
addpath('SuiteSparse/ssget');
addpath('utils');  

tol_pcg = 1e-08;
maxit_pcg = 100;

%% Example using a small matrix with truncations based on full eigendecomposition
% Get a matrix from SuiteSparse
suitesparse_criteria.names = "494_bus";
ids = SuitesSparseHelper.get(suitesparse_criteria);
Prob = ssget(ids(1));  % Prob is a struct (matrix, name, meta-data, ...)
A = Prob.A;
n = size(A, 1);

% incomplete Cholesky
opts.type = 'nofill';
Q = ichol(A, opts);

% Compute preconditioners: here we compare ichol with SVD-based truncation
% and the Bregman log determinant truncation defined in [1]
r = 30;  % rank of low-rank terms

% SVD-based preconditioner
config_svd.method = 'evd';
p_svd = bldp.svd_preconditioner(Q, A, r, config_svd);

% Bregman-based preconditioner
config_breg.method = 'evd';
p_breg = bldp.bregman_preconditioner(Q, A, r, config_breg);

% Solve systems with PCG
b = randn(n, 1);
pcg(A, b, tol_pcg, maxit_pcg);
pcg(A, b, tol_pcg, maxit_pcg, Q, Q');
pcg(A, b, tol_pcg, maxit_pcg, p_svd.action);
pcg(A, b, tol_pcg, maxit_pcg, p_breg.action);

%% Example using a larger matrix using bldp
% Get a matrix from SuiteSparse
suitesparse_criteria.names = "Pres_Poisson";
ids = suitesparse_helper.get(suitesparse_criteria);
Prob = ssget(ids(1));  % Prob is a struct (matrix, name, meta-data, ...)
A = Prob.A;
n = size(A, 1);

% incomplete Cholesky
opts.type = 'nofill';
opts.diagcomp = 1;
Q = ichol(A, opts);

% Compute preconditioners: here we compare ichol with SVD-based truncation
% and the Bregman log determinant truncation defined in [1]
r = 150;  % rank of low-rank terms

% SVD-based preconditioner
config_svd.method = 'nystrom';
p_svd = bldp.svd_preconditioner(Q, A, r, config_svd);

% Bregman-based preconditioner
config_breg.ratio = 0.4;
% We set the `ratio` field to denote how much of the rank `r` budget that
% we allocate towards estimating the largest eigenvalues. That is, we
% estimate `ratio*r` of the largest positive eigenvalues and `(1-ratio)*r`
% of the smallest negative eigenvalues.
config_breg.method = 'krylov_schur';
config_breg.estimate_largest_with_nystrom = 1;
config_breg.tol = 1e-05;
config_breg.maxit = 30;
config_breg.subspace_dim = 100;
p_breg = bldp.bregman_preconditioner(Q, @(x) A*x, r, config_breg);

% Solve systems with PCG
b = randn(n, 1);
pcg(A, b, tol_pcg, maxit_pcg);
pcg(A, b, tol_pcg, maxit_pcg, Q, Q');
pcg(A, b, tol_pcg, maxit_pcg, p_svd.action);
pcg(A, b, tol_pcg, maxit_pcg, p_breg.action);

