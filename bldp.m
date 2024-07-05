classdef bldp 
    methods (Static = true)

        function result = svd_preconditioner(Q, S, r, config)
            n = size(Q, 1);
            if strcmp(config.method, 'nystrom')
                tic;
                if ~isfield(config, 'Omega')
                    config.Omega = randn(n, r);
                end
                [U, D] = bldp.nystrom(Q, S, config.Omega);
                result.U = U;
                result.D = D - speye(size(D, 1));
                result.ctime = toc;
            elseif strcmp(config.method, 'evd')
                [result.U, result.D] = bldp.truncate_matrix(Q \ S / Q' - eye(n), @(x) abs(x), r);
            else
                error("Invalid `method` field in `config`.")
            end
            smw = @ (x) bldp.SMW(result.U, result.D, result.U', x);
            result.action = @ (x) Q' \ (smw(Q \ x));
        end

        function result = bregman_preconditioner(Q, S, r, config)
            tic
            if strcmp(config.method, 'krylov_schur')
                result = bldp.bpc_krylov_schur(Q, S, r, config);
            elseif strcmp(config.method, 'evd')
                result = bldp.bpc_evd(Q, S, r);
                result.diagnostics.flag = 0;
            else
                error("Invalid `method` field in `config`.")
            end
            result.ctime = toc;
            result.action = @ (x) Q' \ (result.action_inner(Q \ x));
        end

        function result = bpc_krylov_schur(Q, S, r, config)
            if ~isa(S, 'function_handle')
                error("S must be a function handle!");
            end
            opts.issym = 1;
            opts.tol = config.tol;
            opts.maxit = config.maxit;
            opts.disp = 0;
            opts.fail = 'drop';
            opts.p = config.subspace_dim;
            if ~isfield(config, 'v')
                opts.v0 = randn(size(Q, 1), 1);
            end
            sigma = 'largestreal';

            n = size(Q, 1);
            g = @ (x) Q \ (S(Q' \ x));

            result.diagnostics = struct;
            if isfield(config, 'r_largest')
                r_largest = config.r_largest;
            elseif isfield(config, 'ratio')
                r_largest = floor(r * config.ratio);
            else
                error("Need to specify `ratio` or `r_largest` in config.")
            end
            r_smallest = r - r_largest;
            result.diagnostics.ks_flag = 0;

            if r_largest > 0
                % Estimate largest eigenvalues
                if config.estimate_largest_with_nystrom
                    Omega = randn(n, r_largest);
                    [U, D] = bldp.nystrom(Q, S, Omega);
                    D = diag(D)';
                    H_max = D(1,1);
                    result.diagnostics.nc = r_largest;
                else
                    [U, H_max, result.diagnostics.ks_flag] = eigs(g, n, r_largest, sigma, opts);
                    D = diag(H_max)';
                    result.diagnostics.nc = size(U, 2);
                end
            else
                % Estimate largest eigenvalue for negative shift
                H_max = eigs(g, n, 1, sigma, opts);
                U = []; D = [];
                result.diagnostics.nc = 1;
            end
            eig_max = H_max(1, 1);

            % Estimate smallest eigenvalues
            if r_smallest > 0
                shg_defl = @ (x) eig_max*x - g(x);
                [V, Dh, flag] = eigs(shg_defl, n, r_smallest, sigma, opts);
                result.diagnostics.ks_flag = result.diagnostics.ks_flag & flag;
                result.diagnostics.nc = result.diagnostics.nc + size(Dh, 2);
                U = [U, V];
                D = [D,  eig_max - diag(Dh)'];
            end

            if min(D) <= 0
                fprintf("[Krylov-Schur] [MATLAB] Negative eigenvalues!");
                result.diagnostics.flag = 1;
                nonsingular_idx = D>0;
                D = D(nonsingular_idx);
                U = U(:, nonsingular_idx);
            end

            result.D = diag(D - 1);
            result.U = U;
            result.r_largest = r_largest;
            result.r_smallest = r_smallest;
            result.action_inner = @ (x) bldp.SMW(result.U, result.D, result.U', x);
        end

        function result = bpc_evd(Q, S, r)
            bregman_curve = @(x) x - log(1 + x);
            [result.U, result.D] = bldp.truncate_matrix(Q \ S / Q' - eye(size(Q, 1)), bregman_curve, r);
            result.action_inner = @ (x) bldp.SMW(result.U, result.D, result.U', x);
        end

        function [U, D] = truncate_matrix(M, f, r)
            M = full(M);
            [U, D] = eig((M + M')/2);
            [~, i] = sort(f(diag(D)));
            i = i(end-r+1:end);
            U = U(:, i);
            D = D(i, i);
        end

        function v = bregman_divergence(X, Y)
            chol(X);
            chol(Y);
            p = X / Y;
            [~, n] = size(X);
            v = trace(p) - log(det(p)) - n;
        end

        function y = SMW(U, C, V, x)
            inner = inv(C) + V*U;
            y = inner \ (V * x);
            y = x - U * y;
        end

        function [U, S] = nystrom(Q, S, Omega)
            if isa(S, 'function_handle')
                Y = Q \ splitapply(S, Q' \ Omega, 1:size(Omega,2));
            else
                Y = Q \ (S * (Q' \ Omega));
            end
            CC = Omega' * Y;
            try
                Cl = chol(CC);
            catch ME
                error("Nyström fails, inner matrix is not PSD.")
            end
            Yhat = Y / Cl;
            [U, S, ~] = svd(Yhat, 'econ');
            S = S * S;
        end
    end
end