classdef bldp 
    methods (Static = true)

        function result = svd_preconditioner(Q, S, config)
            if strcmp(config.method, 'evd')
                tic;
                [U, D, V] = svd(Q \ S / Q' - eye(size(Q, 1)));
                result.U = U(:, 1:config.r);
                result.D = D(1:config.r, 1:config.r);
                result.V = V(:, 1:config.r);
                result.ctime = toc;
            else
                % convert to handle
                if ~isa(S, 'function_handle')
                    S = @(x) S*x;
                end
                % which Nyström approximation?
                mat_action = @(x) Q \ S(Q' \ x) - x;
                if strcmp(config.method, 'indefinite_nystrom')
                    tic;
                    [result.U, result.D, result.V] = bldp.indefinite_nystrom(mat_action, config.sketching_matrix, config.r);
                    result.ctime = toc;
                elseif strcmp(config.method, 'nystrom')
                    if ~isa(S, 'function_handle')
                        S = @(x) S*x;
                    end
                    tic;
                    [result.U, result.D, result.V] = bldp.nystrom(mat_action, config.sketching_matrix);
                    result.ctime = toc;
                elseif strcmp(config.method, 'krylov_schur')
                    opts.issym = 1;
                    opts.tol = config.tol;
                    opts.maxit = config.maxit;
                    opts.disp = 0;
                    opts.fail = 'drop';
                    opts.p = config.subspace_dimension;
                    [result.U, result.D, result.diagnostics.ks_flag] = eigs(mat_action, size(Q, 1), config.r, 'largestabs', opts);
                    result.diagnostics.nc = size(result.D, 2);
                    result.V = result.U;
                    result.ctime = toc;
                else
                    error("Invalid `method` field in `config`.")
                end
            end
            smw = @ (x) bldp.SMW(result.U, result.D, result.V', x);
            result.action = @ (x) Q' \ (smw(Q \ x));
        end

        function result = bregman_preconditioner(Q, S, config)
            tic
            if strcmp(config.method, 'krylov_schur')
                result = bldp.bpc_krylov_schur(Q, S, config);
            elseif strcmp(config.method, 'evd')
                config.reverse = 0;
                result = bldp.bpc_evd(Q, S, config);
                result.diagnostics.flag = 0;
            else
                error("Invalid `method` field in `config`.")
            end
            result.V = result.U;
            result.ctime = toc;
            result.action = @ (x) Q' \ (result.action_inner(Q \ x));
        end

        function result = reverse_bregman_preconditioner(Q, S, config)
            tic
            if strcmp(config.method, 'evd')
                config.reverse = 1;
                result = bldp.bpc_evd(Q, S, config);
                result.diagnostics.flag = 0;
            else
                error("Invalid `method` field in `config`.")
            end
            result.V = result.U;
            result.ctime = toc;
            result.action = @ (x) Q' \ (result.action_inner(Q \ x));
        end

        function result = bpc_krylov_schur(Q, S, config)
            if ~isa(S, 'function_handle')
                error("S must be a function handle!");
            end
            opts.issym = 1;
            opts.tol = config.tol;
            opts.maxit = config.maxit;
            opts.disp = 0;
            opts.fail = 'drop';
            opts.p = config.subspace_dimension;
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
                r_largest = floor(config.r * config.ratio);
            else
                error("Need to specify `ratio` or `r_largest` in config.")
            end
            r_smallest = config.r - r_largest;
            result.diagnostics.ks_flag = 0;

            if r_largest > 0
                % Estimate largest eigenvalues
                if config.estimate_largest_with_nystrom
                    [U, D] = bldp.nystrom(@ (x) Q \ S(Q' \ x), config.sketching_matrix);
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
                if config.estimate_largest_with_nystrom
                    [~, H_max] = bldp.nystrom(@ (x) Q \ S(Q' \ x), config.sketching_matrix);
                    result.diagnostics.nc = 1;
                else
                    H_max = eigs(g, n, 1, sigma, opts);
                end
                U = []; D = [];
                result.diagnostics.nc = 1;
            end
            if isempty(H_max)
                error("Failed to estimate largest eigenvalues. " + ...
                      "Increase max iterations or subspace dimension.");
            end

            % Estimate smallest eigenvalues
            if r_smallest > 0
                eig_max = H_max(1, 1);
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

        function result = bpc_evd(Q, S, config)
            if config.reverse
                f = @(x) 1./(1+x) + log(1 + x) - 1;
            else
                f = @(x) x - log(1 + x);
            end
            [result.U, result.D] = bldp.truncate_matrix(Q \ S / Q' - eye(size(Q, 1)), f, config.r);
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
            p = X / Y;
            v = trace(p) - log(det(p)) - size(X, 1);
        end

        function y = SMW(U, C, V, x)
            inner = inv(C) + V*U;
            y = inner \ (V * x);
            y = x - U * y;
        end

        function [U, S, V] = nystrom(A, sketching_matrix)
            if ~isa(A, 'function_handle')
                error("A must be a function handle.");
            end
            Y = splitapply(A, sketching_matrix, 1:size(sketching_matrix, 2));
            inner = sketching_matrix' * Y;
            try
                C = chol(inner);
                Yhat = Y / C;
                [U, S] = svd(Yhat, 'econ');
                S = S * S;
                V = U;
            catch
                [U, S, V] = bldp.indefinite_nystrom(A, sketching_matrix, size(sketching_matrix, 2));
            end
        end

        function [U, S, V] = indefinite_nystrom(A, sketching_matrix, r)
            % Implements:
            % Nakatsukasa, Yuji, and Taejun Park. "Randomized low-rank
            % approximation for symmetric indefinite matrices." SIAM 
            % Journal on Matrix Analysis and Applications 44.3 (2023):
            % 1370-1392.
            if ~isa(A, 'function_handle')
                error("A must be a function handle.");
            end
            Y = splitapply(A, sketching_matrix, 1:size(sketching_matrix, 2));
            [U, S, V] = svd(sketching_matrix' * Y);
            U = Y * U(:, 1:r);
            S = pinv(S(1:r, 1:r));
            V = Y * V(:, 1:r);
        end
    end
end