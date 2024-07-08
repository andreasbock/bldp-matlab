clear; clearvars; close all; beep off;
rng("default")

% Generate a random symmetric matrix
n = 150;
A = randn(n, n);

%% Positive definite matrix
A = A*A';
A_hdl = @(x) A*x;

% Sanity check: no low-rank approximation
r = n;

% Sample Gaussian
Omega = randn(n, r);

% Nyström (naive)
Y = A * Omega;
A_nys_naive = Y * ((Omega' * Y) \ Y');

% Nyström (bldp)
[U, S] = bldp.nystrom(A_hdl, Omega);
A_nys_bldp = U * S * U';

% Indefinite Nyström (c = 1)
[U, S, V] = bldp.indefinite_nystrom(A_hdl, Omega, r);
A_nys_bldp_indef = U * S * V';

% Random SVD (single pass)
[O, ~] = qr(Y);
O = O(:, 1:r);
Pi = (O' * Y) / (O' * Omega);
A_rsvd_sp = O*Pi*O';

error_sanity_check = max([ ...
    norm(A - A_nys_naive) ...
    norm(A - A_nys_bldp) ...
    norm(A - A_nys_bldp_indef) ...
    norm(A - A_rsvd_sp)] ...
);
if error_sanity_check > 1e-08
    error("PSD matrix: error in Nyström approximations!");
end

% r = n * 0.5
r = ceil(n/ 2);

% Sample Gaussian
Omega = randn(n, r);
Y = A * Omega;

% Nyström (naive)
A_nys_naive = Y * ((Omega' * Y) \ Y');

% Nyström (bldp)
[U, S] = bldp.nystrom(A_hdl, Omega);
A_nys_bldp = U * S * U';

% Indefinite Nyström (c = 1)
[U, S, V] = bldp.indefinite_nystrom(A_hdl, Omega, r);
A_nys_bldp_indef = U * S * V';

err_nys_bldp = norm(A_nys_naive - A_nys_bldp);
err_nys_bldp_indef = norm(A_nys_naive - A_nys_bldp_indef);
err_nys_impl = max(err_nys_bldp, err_nys_bldp_indef);

if err_nys_impl > 1e-08
    error("bldp implementation differs from Nystrom!");
end

eA = real(eig(A));
eA_rsvd_sp = real(eig(A_rsvd_sp));
eA_nys = real(eig(A_nys_bldp_indef));

max(eA), min(eA)
max(eA_nys), min(eA_nys)


%% INDEFINITE MATRIX
[V, D] = eig(A);
flip_sign = (rand(n, 1)<.5)*2 - 1;
D = diag(diag(D) .* flip_sign);
A = V * D * V';
A_hdl = @(x) A*x;

% Sanity check: no low-rank approximation
r = n;

% Sample Gaussian
Omega = randn(n, r);

% Nyström (naive)
Y = A * Omega;
A_nys_naive = Y * ((Omega' * Y) \ Y');

% Nyström (bldp)
[U, S, V] = bldp.nystrom(A_hdl, Omega);
A_nys_bldp = U * S * V';

% Indefinite Nyström (c = 1)
[U, S, V] = bldp.indefinite_nystrom(A_hdl, Omega, r);
A_nys_bldp_indef = U * S * V';

% Random SVD (single pass)
[O, ~] = qr(Y);
O = O(:, 1:r);
Pi = (O' * Y) / (O' * Omega);
A_rsvd_sp = O*Pi*O';

error_sanity_check = max([ ...
    norm(A - A_nys_naive) ...
    norm(A - A_nys_bldp) ...
    norm(A - A_nys_bldp_indef) ...
    norm(A - A_rsvd_sp)] ...
);
if error_sanity_check > 1e-08
    error("Indefinite matrix: error in Nyström approximations!");
end

% r = n * 0.5
r = ceil(n/ 2);

% Sample Gaussian
c = 1.5;
cr = round(c*r);
Omega = randn(n, cr);
Y = A * Omega;
Omega_r = Omega(:, 1:r);

% Truncated SVD
[A_U, A_L, A_V] = svd(A, 'econ');
Ar = A_U(:,1:r)*A_L(1:r,1:r)*A_V(:,1:r)';

% Indefinite Nyström (c = 1.5)
[U, S, V] = bldp.indefinite_nystrom(A_hdl, Omega, r);
A_nys_bldp_indef = U * S * V';

% Random SVD (single pass)
[O, ~] = qr(Y);
O = O(:, 1:r);
Pi = (O' * Y) / (O' * Omega);
A_rsvd_sp = O*Pi*O';

err_nys_bldp_indef = norm(Ar - A_nys_bldp_indef);
err_rsvd = norm(Ar - A_rsvd_sp);

eA = real(eig(A));
eAr = real(eig(Ar));
eA_rsvd_sp = real(eig(A_rsvd_sp));
eA_nys = real(eig(A_nys_bldp_indef));

max(eA), min(eA)
max(eA_nys), min(eA_nys)


plot(1:n, sort(eA_nys)); hold on;
plot(1:n, sort(eA_rsvd_sp), 'Color', 'red'); hold on;
plot(1:n, sort(Ar), 'Color', 'black'); hold off;


% TODO: finish testing somehow