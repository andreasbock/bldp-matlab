clear; clearvars; close all; beep off;

% for suitesparse_helper
addpath('utils');  

tol_pcg = 1e-08;
maxit_pcg = 200;

% Example using a small matrix with truncations based on full eigendecomposition
% Fetch a matrix from SuiteSparse
suitesparse_criteria.names = "494_bus";
ids = SuitesSparseHelper.get(suitesparse_criteria);
Prob = ssget(ids(1));  % Prob is a struct (matrix, name, meta-data, ...)
A = Prob.A;

% incomplete Cholesky
opts.type = 'nofill';
Q = ichol(A, opts);

% Compute preconditioners: here we compare ichol with SVD-based truncation
% and the Bregman log determinant truncation defined in [1]
r = 30;  % rank of low-rank terms

config.method = 'evd';
config.r = r;

% SVD-based preconditioner
p_svd = bldp.svd_preconditioner(Q, A, config);

% Bregman-based preconditioner
p_breg = bldp.bregman_preconditioner(Q, A, config);

% Solve systems with PCG
b = randn(size(A, 1), 1);
pcg(A, b, tol_pcg, maxit_pcg);
pcg(A, b, tol_pcg, maxit_pcg, Q, Q');
pcg(A, b, tol_pcg, maxit_pcg, p_svd.action);
pcg(A, b, tol_pcg, maxit_pcg, p_breg.action);

%% Example using a larger matrix using bldp
% Get a matrix from SuiteSparse
suitesparse_criteria.names = "Pres_Poisson";
%suitesparse_criteria.names = "crankseg_1";

ids = SuitesSparseHelper.get(suitesparse_criteria);
Prob = ssget(ids(1));  % Prob is a struct (matrix, name, meta-data, ...)
A = Prob.A;
n = size(A, 1);

% incomplete Cholesky
opts.type = 'ict';
opts.diagcomp = 0;
opts.droptol = 1e-03;
Q = ichol(A, opts);

% Compute preconditioners: here we compare ichol with SVD-based truncation
% and the Bregman log determinant truncation defined in [1]
r = 325;  % rank of low-rank terms

% Park-Nakatsukasa Nyström (i.e. Nyström for indefinite matrices)
c = 1.5;  % oversampling as a percentage of rank
config_nys_indef.r = r;
config_nys_indef.sketching_matrix = randn(n, round(r*c));
config_nys_indef.method = 'indefinite_nystrom';
p_nys_indef = bldp.svd_preconditioner(Q, A, config_nys_indef);

% Bregman-based preconditioner
config_breg.r = r;
config_breg.ratio = 0.8;
% We set the `ratio` field to denote how much of the rank `r` budget that
% we allocate towards estimating the largest eigenvalues. That is, we
% estimate `ratio*r` of the largest positive eigenvalues and `(1-ratio)*r`
% of the smallest negative eigenvalues.
config_breg.method = 'krylov_schur';
config_breg.estimate_largest_with_nystrom = 0;
config_breg.tol = 1e-01;
config_breg.maxit = 5;
config_breg.subspace_dimension = r + 5;  % must be larger than `r`
p_breg = bldp.bregman_preconditioner(Q, @(x) A*x, config_breg);

% Solve systems with PCG
b = randn(n, 1);
pcg(A, b, tol_pcg, maxit_pcg);
pcg(A, b, tol_pcg, maxit_pcg, Q, Q');
pcg(A, b, tol_pcg, maxit_pcg, p_nys_indef.action);
pcg(A, b, tol_pcg, maxit_pcg, p_breg.action);

