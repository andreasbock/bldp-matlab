%% Test Sherman-Morrison-Woodbury identity
n = 10;
r = 5;

Q = randn(n, r);
O = qr(Q);
S = diag(abs(randn(1, r)));
A = O*S*O';

in = randn(n, 1);
out1 = (eye(n) + A) \ in;
out2 = bldp.SMW(O, S, O', in);

err = norm(out1 - out2);
if err > 1e-14
    error("Test failed");
else
    disp("Test passed.")
end