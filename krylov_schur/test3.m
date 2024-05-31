clear; clearvars; close all; beep off; format long;
global A;

N = 30;       % size of matrix
k = 5;        % number of eigenvalues needed
m = 10;       % Krylov restart threshold
maxIt = 20;
tol = 1e-05;

B = rand(N, N); B = B*B';
[ev_B, e_B] = eig(B);
e_B = diag(e_B);

eps = max(e_B) + 100;
A = eps*eye(N) - B;
[ev, es] = eig(A);

es = diag(es);
[~, ix] = sort(abs(es), 'descend');
es = es(ix);
disp('eigenvalues by matlab command eig() '); disp(es(1:k));

v1 = rand(N, 1);
[Q, H, isC, flag_, nc_, ni_] = KrylovSchur(@Ax, v1, N, k, m, maxIt, tol);
[e_ks, v, ~, flag, nc, ni] = KrylovSchurEig(@Ax, v1, N, k, m, maxIt, tol);

if nc ~= k
    error('KS did not converge.');
end
[~, ix] = sort(abs(e_ks), 'descend');
e_ks = e_ks(ix);

disp('eigenvalues '); disp(e_B);
e_ks_epsed = eps - e_ks;
disp('eigenvalues by KrylovSchurEig '); disp(e_ks_epsed);

disp(['error of KrylovSchurEig : ', num2str(norm(A * v - v * diag(e_ks)))]);
fprintf(1, 'number of desired eigenvalues: %d \n', k);
fprintf(1, 'number of converged eigenvalues: %d \n', nc);
fprintf(1, 'number of iteration used : %d \n\n', ni);

sz = 50;
scatter((e_ks_epsed), 1.2*ones(1, length(e_ks_epsed)), sz, 'black', 'filled'); hold on;
scatter((e_B), ones(1, length(e_B)), sz, 'red', 'filled'); hold on;
pbaspect([10 1 1])
hold off;
