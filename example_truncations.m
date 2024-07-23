% This script constructs the Bregman log determinant truncation and
% compares it to an SVD in the context of preconditioning a symmetric
% linear system of equations.

% For ease of exposition some computational inefficiencies are present.
% See `example_bldp.m` where we show similar examples using the proposed
% `bldp.m` library with more scalable/realistic implementations of what
% goes on in this script.

addpath('SuiteSparse/ssget');
addpath('utils');

% Set up Sx=b
% Get a matrix from SuiteSparse
suitesparse_criteria.names = "494_bus";
ids = SuitesSparseHelper.get(suitesparse_criteria);
S = ssget(ids(1)).A;

% Compute G such that S == Q(I + G)Q'
n = size(S, 1);
I = eye(n);
Q = ichol(S);
G = Q\S/Q' - I;
G = (G + G')/2;

% set rank
r = 30;

% SVD of G
[V, E, W] = svd(G);
V_svd = V(:, 1:r);
W_svd = W(:, 1:r);
E_svd = E(1:r, 1:r);
G_svd = V_svd * E_svd * W_svd';

% Bregman log determinant truncation of G
f = @(x) x - log(1 + x);
[U, D] = eig(G);
[~, i] = sort(f(diag(D)));
idx_bld = i(end-r+1:end);
U_bld = U(:, idx_bld);
D_bld = D(idx_bld, idx_bld);
G_bld = U_bld * D_bld * U_bld';

% Construct SVD preconditioner
% (bldp.SMW just implements the Sherman-Morrison-Woodbury identity)
inner_svd = @(x) bldp.SMW(V_svd, E_svd, W_svd', x);
P_svd = @(x) Q' \ (inner_svd(Q \ x));

% Construct Bregman log determinant preconditioner
inner_bld = @(x) bldp.SMW(U_bld, D_bld, U_bld', x);
P_bld = @(x) Q' \ (inner_bld(Q \ x));

% Solve systems using PCG
b = randn(n, 1);
tol = 1e-08;
maxit = 200;
pcg(S, b, tol, maxit);
pcg(S, b, tol, maxit, Q, Q');
pcg(S, b, tol, maxit, P_svd);
pcg(S, b, tol, maxit, P_bld);

