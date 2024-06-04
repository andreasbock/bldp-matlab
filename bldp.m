classdef bldp 
    methods (Static = true)
        
        function parse_config(config)
            if ~(config.nystrom + config.krylov_schur + config.evd) == 1
                error("Config is invalid. Must specify one of" + ...
                    " config.nystrom, config.krylov_schur or config.evd!")
            end
        end

        function result = svd_preconditioner(Q, S, r, config)
            arguments
                Q; S; r; config = struct;
            end
            bldp.parse_config(config);

            n = size(Q, 1);
            if config.nystrom
                % for large n we use a stabilised Nyström:
                rc = round(r * config.c);
                % Nakatsukasa, Yuji, and Taejun Park. "Randomized low-rank approximation for symmetric
                % indefinite matrices." SIAM Journal on Matrix Analysis and Applications 44.3 (2023): 1370-1392.
                Omega = randn(n, rc);
                Y = Q \ (S * (Q' \ Omega)) - Omega;
                inner = Omega' * Y;
                [V_inner, D_inner] = eig((inner + inner') / 2);
                [~, i] = sort(abs(real(diag(D_inner))));
                i = i(end-r+1:end);
                iinner_truncated = V_inner(:, i) * pinv(D_inner(i, i)) * V_inner(:, i)';
                if config.full_assembly
                    nys = Y * iinner_truncated * Y';
                    result.nys = nys;
                    [V, D] = eig((nys + nys')/2);
                    [~, i] = sort(abs(real(diag(D))));
                    i = i(end-r+1:end);
                    result.V = V(:, i);
                    result.D = D(i, i);
                    if min(diag(result.D)) <= -1
                        error("Nyström preconditioner is singular.")
                    end
                end
                smw = @ (x) bldp.SMW(Y, iinner_truncated, Y', x);
            elseif config.krylov_schur
                addpath('krylov_schur');
                restart = max(config.restart, r) + 1;
                [O, H, ~, flag, ~, ni] = KrylovSchur(@ (x) Q \ (S * (Q' \ x)) - x, config.v, n, r, restart, config.maxit, config.tol);
                result.V = O(:, 1:r);
                result.D = diag(diag(H));
                smw = @ (x) bldp.SMW(U, L, U', x);
                result.diagnostics.ni = ni;
                result.diagnostics.flag = flag;
            elseif config.evd
                G = full(Q \ S / Q' - eye(n));
                [V, D] = eig((G + G')/2);
                [~, i] = sort(abs(real(diag(D))));
                i = i(end-r+1:end);
                result.V = V(:, i);
                result.D = D(i, i);
                smw = @ (x) bldp.SMW(result.V, result.D, result.V', x);
            else
                error("Invalid config.")
            end
            result.action = @ (x) Q' \ (smw(Q \ x));
        end

        function result = bregman_preconditioner(Q, S, r, config)
            arguments
                Q; S; r; config = struct;
            end
            result = bldp.bpc(Q, S, r, @(x) x - 1 - log(x), config);
        end

        function result = reverse_bregman_preconditioner(Q, S, r, config)
             arguments
                 Q; S; r; config = struct;
             end
             result = bldp.bpc(Q, S, r, @(x) 1./(1 + x) - log(1./(1 + x)) - 1, config);
        end

        function result = bpc(Q, S, r, f, config)
            if ~config.krylov_schur
                result = bldp.bpc_evd(Q, S, r, f);
                result.diagnostics.flag = 0;
            else
                result = bldp.bpc_krylov_schur(Q, S, r, f, config);
            end
            result.action = @ (x) Q' \ (result.action_inner(Q \ x));
        end

        function result = bpc_krylov_schur(Q, S, r, ~, config)
            addpath('krylov_schur');
            n = size(S, 1);
            g = @ (x) Q \ (S * (Q' \ x));

            r_largest = ceil(r * config.bregman.ratio);
            r_smallest = r - r_largest;
            restart = max(config.restart, r) + 1;

            if r_largest > 0
                % Estimate largest eigenvalues
                [O, H_max, ~, diagnostics.flag, ~, diagnostics.ni] = KrylovSchur(g, config.v, n, r_largest, restart, config.maxit, config.tol);
                U = [O(:, 1:r_largest)];
                D = diag(H_max(1:r_largest, 1:r_largest))';
            else
                % Estimate largest eigenvalue for negative shift
                [O, H_max, ~, diagnostics.flag, ~, diagnostics.ni] = KrylovSchur(g, config.v, n, 1, restart, config.maxit, config.tol);
                U = []; D = [];
            end
            eig_max = H_max(1, 1);

            % Estimate smallest eigenvalues
            if r_smallest > 0
                shg_defl = @ (x) eig_max*x - g(x);
                [O, H_min, ~, flag, ~, ni] = KrylovSchur(shg_defl, bldp.lincomb(O), n, r_smallest, restart, config.maxit, config.tol);
                
                U = [U, O(:, 1:r_smallest)];
                D = [D,  eig_max - diag(H_min(1:r_smallest, 1:r_smallest))'];

                diagnostics.ni = diagnostics.ni + ni;
                diagnostics.flag = diagnostics.flag || flag;
            end

            result.D = diag(D - 1);
            result.U = U;
            result.r_largest = r_largest;
            result.r_smallest = r_smallest;
            result.action_inner = @ (x) bldp.SMW(result.U, result.D, result.U', x);
            result.diagnostics = diagnostics;
        end

        function result = bpc_evd(Q, S, r, f)
            G = Q \ S / Q';
            G = (G + G') / 2;
            [U, E] = eig(full(G));
            eig_G = real(diag(E));
            [~, idx] = sort(f(eig_G));
            idx_r = idx(end-r+1:end);
            if norm(imag(f(eig_G))) > 0 || norm(imag(eig_G)) > 0
                error('Imagininary eigenvalues in approximations.')
            end
            result.idx = idx_r;
            result.U = U(:, idx_r);
            result.D = diag(eig_G(idx_r) - 1);
            result.action_inner = @ (x) bldp.SMW(result.U, result.D, result.U', x);
        end

        function v = bregman_divergence(X, Y)
            chol(X);
            chol(Y);
            p = X / Y;
            [~, n] = size(X);
            v = trace(p) - log(det(p)) - n;
        end

        function [mu1, mu2] = extremal_eigenvalues(A, subspace_dimension)
                e = eigs(A, 2, 'bothendsreal', 'SubspaceDimension', subspace_dimension);
                mu1 = e(1);
                mu2 = e(2);
        end

        function y = SMW(U, C, V, x)
            inner = inv(C) + V*U;
            y = inner \ (V * x);
            y = x - U * y;
        end

        function v = lincomb(V)
            n = size(V, 2);
            v = sum(V, 2) ./ n;
        end
    end
end