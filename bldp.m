classdef bldp 
    methods (Static = true)
        
        function n = large_matrix_threshold()
            n = 5001;
        end

        function [mu1, mu2] = extremal_eigenvalues(A, subspace_dimension)
                e = eigs(A, 2, 'bothendsreal', 'SubspaceDimension', subspace_dimension);
                mu1 = e(1);
                mu2 = e(2);
        end

        function y = ShermanMorrisonWoodburyIdentity(U, C, V, x)
            inner = inv(C) + V*U;
            y = inner \ (V * x);
            y = x - U * y;
        end

        function [PC_action, PC, G_lowrank, eigenvalues, nearness] = svd_preconditioner(Q, S, r)
            n = size(Q, 1);
            if n < bldp.large_matrix_threshold()
                I = eye(n);
                G = full(Q \ S / Q' - I);
                [U, Lambda, V] = svd(G);
                G_lowrank = U(:,1:r) * Lambda(1:r,1:r) * V(:,1:r)';
                smw = @ (x) bldp.ShermanMorrisonWoodburyIdentity(U(:,1:r), Lambda(1:r,1:r), V(:,1:r)', x);
                PC = Q * (eye(n) + G_lowrank) * Q';
                eigenvalues = eig(G);
                nearness.cond = condest((I + G_lowrank)\(I + G));
                nearness.div = bldp.bregman_divergence(I + G, I + G_lowrank);
            else
                % for large n we use Nyström
                Y = Omega + Q \ (S * (Q' \ Omega));
                smw = @ (x) bldp.ShermanMorrisonWoodburyIdentity(Y, inv(Omega'*Y), Y', x);
                eigenvalues = 0;
                PC = 0;
                G_lowrank = 0;
                nearness.cond = 0;
                nearness.div = 0;
            end
            PC_action = @ (x) Q' \ (smw(Q \ x));
        end

        function [PC_action, PC, G_lowrank, eigenvalues, nearness] = bregman_preconditioner(Q, S, r)
             %f = @(x) x - log(1 + x);
             f = @(x) x - 1 - log(x);
             [PC_action, PC, G_lowrank, eigenvalues, nearness] = bldp.bpc(Q, S, r, f);
        end

        function [PC_action, PC, G_lowrank, eigenvalues, nearness] = reverse_bregman_preconditioner(Q, S, r)
             f = @(x) 1./(1 + x) - log(1./(1 + x)) - 1;
             [PC_action, PC, G_lowrank, eigenvalues, nearness] = bldp.bpc(Q, S, r, f);
        end

        function [d, V] = sorted_eig(M)
            [W, E] = eig(full(M));
            [d, ind] = sort(diag(E), 'descend');
            d = real(d);
            V = real(W(:,ind));
        end

        function [PC, G_compensated, G, V, L] = bregman_preconditioner2(Q, S, r)
            [n, ~] = size(Q);
            G_no1 = Q \ S / Q';
            f = @(x) x - 1 - log(x);
            [lambda, Basis] = bldp.sorted_eig(G_no1);
            [~, idx] = sort(f(lambda));
            idx_r = idx(end-r+1:end);
            V = Basis(:, idx_r);
            lambda = double(lambda(idx_r));
            G = G_no1 - eye(n);
            lambda = lambda - double(1);

            L = diag(lambda);
            G_compensated = V * L * V';
            PC = Q * (eye(n) + G_compensated) * Q';
        end

        function [PC_action, PC, G_lowrank, eigenvalues, nearness] = bpc(Q, S, r, f)
            % Provides a Bregman preconditioner (sorting eigenvalues by `f`).
            n = size(Q, 1);
            if n < bldp.large_matrix_threshold()
                I = eye(n);
                G = full(Q \ S / Q');
                [eig_G, U] = bldp.sorted_eig(G);

                [~, idx] = sort(f(eig_G));
                idx_r = idx(end-r+1:end);
                if norm(imag(f(eig_G))) > 0 || norm(imag(eig_G)) > 0
                    error('Imagininary eigenvalues in approximations.')
                end

                eigenvalues = eig_G(idx_r) - 1;
                U = U(:, idx_r);
                L = diag(eigenvalues);
                smw = @ (x) bldp.ShermanMorrisonWoodburyIdentity(U, L, U', x);
                G_lowrank = U * L * U';
                PC = Q * (eye(n) + G_lowrank) * Q';
                
                nearness.cond = condest((I + G_lowrank)\(I + G));
                nearness.div = bldp.bregman_divergence(I + G, I + G_lowrank);
            else
                nearness.cond = 0;
                nearness.div = 0;
                error("TODO!")
            end
            PC_action = @ (x) Q' \ (smw(Q \ x));
        end

        function v = bregman_divergence(X, Y)
            chol(X);
            chol(Y);
            p = X / Y;
            [~, n] = size(X);
            v = trace(p) - log(det(p)) - n;
        end

        function plot_bregman_curves(e, e_r, e_apx, path, log_scale, legends)
            font_size = 16;
            mkdir(path);
            if nargin <= 6
                legends = {"$S$", "$\hat{S}_r$ (TSVD)", "$\bar{S}_r$ (BT)", "Overlap"};%, '', '', ''];
            end
            legends = {'$G$', '$\left[ \! \left[ G\right] \! \right]_r$ (TSVD)', '$\langle\!\langle G\rangle\!\rangle_r$ (BT)', 'Overlap'};
            config = Config();
            svd_colour = config.svd.colour;
            breg_colour = config.breg.colour;
            both_colour = '#379c37';
            alpha_og = 1;
            alpha_r = 1;
            alpha_breg = 1;
            alpha_both = 1;

            % plot Bregman curve
            sz_breg = 150;
            sz_svd = 150;
            sz_both = 140;
            sz_original = 50; 
            figure('Visible', 'off');
            t = @ (x) x - log(1 + x);

            e_both = [];
            e_r_filtered = [];
            e_apx_filtered = [];
            e_filtered = [];
            r = numel(e_r);
            for i = 1:r
                if bldp.flismember(e_r(i), e_apx)
                    e_both = [e_both e_r(i)];
                else
                    e_r_filtered = [e_r_filtered e_r(i)];
                end
            end
            for i = 1:r
                if ~ bldp.flismember(e_apx(i), e_both)
                    e_apx_filtered = [e_apx_filtered e_apx(i)];
                end
            end

            for i = 1:numel(e)
                bool_both = bldp.flismember(e(i), e_both);
                bool_r = bldp.flismember(e(i), e_r_filtered);
                bool_apx = bldp.flismember(e(i), e_apx_filtered);
                if ~ (bool_both || bool_r || bool_apx)
                    e_filtered = [e_filtered e(i)];
                end
            end
            if numel(e_r_filtered) ~= numel(e_apx_filtered)
                error('We have numel(e_r_filtered) ~= numel(e_apx_filtered)');
            end
            
            if numel(e_r_filtered) > 0
                yline(t(1 + e_r_filtered), 'LineStyle', '-', 'LineWidth', 1, 'Alpha', 1, 'Color', svd_colour); hold on;
            end
            if numel(e_apx_filtered) > 0
                yline(t(1 + e_apx_filtered), 'LineStyle', ':', 'LineWidth', 2, 'Alpha', 1, 'Color', breg_colour); hold on;
            end
            
            %scatter(e, zeros(1, numel(e)), sz_original, 'black', 'filled', 'o', 'MarkerEdgeAlpha', 0.0, 'MarkerFaceAlpha', 0.8); hold on;
            s1 = scatter(e_apx_filtered, t(e_apx_filtered), sz_breg, 'filled', 'o', 'MarkerEdgeAlpha', .10, 'MarkerFaceAlpha', alpha_breg, 'MarkerEdgeColor', breg_colour, 'MarkerFaceColor', breg_colour); hold on;
            s2 = scatter(e_r_filtered, t(e_r_filtered), sz_svd, 'filled', 'd', 'MarkerEdgeAlpha', 1.0, 'MarkerFaceAlpha', alpha_r, 'MarkerEdgeColor', '#00ecec', 'MarkerFaceColor', svd_colour); hold on;
            s3 = scatter(e_filtered, t(e_filtered), sz_original, 'black', 'filled', 'o', 'MarkerEdgeAlpha', 0.0, 'MarkerFaceAlpha', alpha_og); hold on;
            s4 = scatter(e_both, t(e_both), sz_both, 'filled', 'o', 'MarkerEdgeAlpha', 0.0, 'MarkerFaceAlpha', alpha_both, 'MarkerEdgeColor', both_colour, 'MarkerFaceColor', both_colour); hold on;
            
            if log_scale
                set(gca, 'YScale', 'log');
            end
            if 0
                for i = 1:numel(e)
                    ei = e(i);
                    ylo = t(ei); 
                    if ylog
                        ylo = log(ylo); 
                    end
                    plot([ei ei], [ylo 0.1],'DisplayName', '', 'LineWidth', 0.4, 'Color', 'black', 'LineStyle', ':')
                end
            end
            xs = linspace(min(e), max(e)+1e-01, 1000);
            plot(xs, t(xs), 'black'); hold on;
            xlabel('$\lambda$', 'Interpreter', 'latex', 'FontSize', font_size)
            %ylabel('$\frac{1}{1+\lambda}-\log(\frac{1}{1+\lambda})-1$', 'Interpreter', 'latex', 'FontSize', font_size)
            ylabel('$\lambda-\log(1+\lambda)$', 'Interpreter', 'latex', 'FontSize', font_size)
            axis square;
            hold off;
            ldg = legend([s3 s2 s1 s4], legends);
            set(ldg, 'Interpreter', 'latex', 'FontSize', font_size, 'Location', 'northoutside', 'Orientation', 'horizontal', 'NumColumns', 4);
            grid off;
            filename = 'bregman_curve.pdf';
            exportgraphics(gcf, fullfile(path, filename));
        end
      
        function plot_svd_curve(e, e_r, e_apx, path, ylog)
            font_size = 16;
            config = Config();
            svd_colour = config.svd.colour;
            breg_colour = config.breg.colour;
            both_colour = '#379c37';
            alpha_og = 1;
            alpha_r = 1;
            alpha_breg = 1;
            alpha_both = 1;

            % plot Bregman curve
            sz_breg = 150;
            sz_svd = 150;
            sz_both = 140;
            sz_original = 100;
            shf = 0;
            figure('Visible', 'off');
            t = @(x) abs(x);
            e_both = [];
            e_r_filtered = [];
            e_apx_filtered = [];
            e_filtered = [];
            r = numel(e_r);
            for i = 1:r
                if bldp.flismember(e_r(i), e_apx)
                    e_both = [e_both e_r(i)];
                else
                    e_r_filtered = [e_r_filtered e_r(i)];
                end
            end
            for i = 1:r
                if ~ bldp.flismember(e_apx(i), e_both)
                    e_apx_filtered = [e_apx_filtered e_apx(i)];
                end
            end

            for i = 1:numel(e)
                bool_both = bldp.flismember(e(i), e_both);
                bool_r = bldp.flismember(e(i), e_r_filtered);
                bool_apx = bldp.flismember(e(i), e_apx_filtered);
                if ~ (bool_both || bool_r || bool_apx)
                    e_filtered = [e_filtered e(i)];
                end
            end
            if numel(e_both) == numel(e_r)
                2;%error('e_tsvd == e_bt');
            end
            if numel(e_r_filtered) ~= numel(e_apx_filtered)
                error('numel(e_r_filtered) ~= numel(e_apx_filtered)');
            end
            if numel(e_r_filtered) ~= 0
                yline(t(shf + e_r_filtered), 'LineStyle', '-', 'LineWidth', 1, 'Alpha', 0.3, 'Color', svd_colour); hold on;
            end
            if numel(e_apx_filtered) ~= 0
                yline(t(shf + e_apx_filtered), 'LineStyle', ':', 'LineWidth', 2, 'Alpha', 0.1, 'Color', breg_colour); hold on;
            end
            
            xs = linspace(min(e), max(e)+1e-01, 100); hold on;
            plot(xs, arrayfun(t, xs), 'black');
            
            e_zeros = zeros(1, numel(e));
            %scatter(e, e_zeros, sz_original, 'black', 'filled', 'o', 'MarkerEdgeAlpha', 0.0, 'MarkerFaceAlpha', 0.8); hold on;
            s1 = scatter(shf + e_apx_filtered, t(shf + e_apx_filtered), sz_breg, 'filled', 'o', 'MarkerEdgeAlpha', 0.0, 'MarkerFaceAlpha', alpha_breg, 'MarkerEdgeColor', breg_colour, 'MarkerFaceColor', breg_colour); hold on;
            s2 = scatter(shf + e_r_filtered, t(shf + e_r_filtered), sz_svd, 'filled', 'd', 'MarkerEdgeAlpha', 1.0, 'MarkerFaceAlpha', alpha_r, 'MarkerEdgeColor', '#00ecec', 'MarkerFaceColor', svd_colour); hold on;
            s3 = scatter(shf + e_filtered, t(shf + e_filtered), sz_original, 'black', 'filled', 'o', 'MarkerEdgeAlpha', 0.0, 'MarkerFaceAlpha', alpha_og);
            s4 = scatter(shf + e_both, t(shf + e_both), sz_both, 'filled', 'o', 'MarkerEdgeAlpha', 0.0, 'MarkerFaceAlpha', alpha_both, 'MarkerEdgeColor', both_colour, 'MarkerFaceColor', both_colour); hold on;
            
            if ylog
                set(gca, 'YScale', 'log');
            end
            if 0
                for i = 1:numel(e)
                    ei = e(i);
                    ylo = t(ei); 
                    if ylog
                        ylo = log(ylo); 
                    end
                    plot([ei ei], [ylo 0.1],'DisplayName', '', 'LineWidth', 0.4, 'Color', 'black', 'LineStyle', ':')
                end
            end
            pause(0.1);
            xlabel('$\lambda$', 'Interpreter', 'latex', 'FontSize', font_size)
            ylabel('$|\lambda|$', 'Interpreter', 'latex', 'FontSize', font_size)
            axis square;
            hold off;
            ldg = legend([s3 s2 s1 s4], {'$G$', '$\left[ \! \left[ G\right] \! \right]_r$ (TSVD)', '$\langle\!\langle G\rangle\!\rangle_r$ (BT)', 'Overlap'});
            set(ldg, 'Interpreter', 'latex', 'FontSize', font_size, 'Location', 'northoutside', 'Orientation', 'horizontal', 'NumColumns', 4);%, 'Location', 'northoutside', 'Orientation', 'horizontal', 'NumColumns', 4);
            grid off;
            filename = 'svd_curve.pdf';
            exportgraphics(gcf, fullfile(path, filename));
        end
       
        function b = flismember(x, xs)
            b = 0;
            for i = 1:numel(xs)
               if abs(x - xs(i)) < 1e-13
                   b = 1;
               end
            end
        end

    end
end