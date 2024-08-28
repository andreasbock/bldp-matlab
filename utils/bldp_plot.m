classdef bldp_plot 
    methods (Static = true)
        
        function plot_bregman_curves(e, e_r, e_apx, path, log_scale, legends)
            font_size = 16;
            mkdir(path);
            %legends = {'$G$', '$\left[ \! \left[ G\right] \! \right]_r$ (TSVD)', '$\langle\!\langle G\rangle\!\rangle_r$ (BLD)', 'Overlap'};
            legends = {'$\left[ \! \left[ G\right] \! \right]_r$ (TSVD)', '$\langle\!\langle G\rangle\!\rangle_r$ (BLD)', 'Overlap'};
            config = Plotting();
            svd_colour = config.svd.colour;
            breg_colour = config.breg.colour;

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
                if bldp_plot.flismember(e_r(i), e_apx)
                    e_both = [e_both e_r(i)];
                else
                    e_r_filtered = [e_r_filtered e_r(i)];
                end
            end
            for i = 1:r
                if ~ bldp_plot.flismember(e_apx(i), e_both)
                    e_apx_filtered = [e_apx_filtered e_apx(i)];
                end
            end

            for i = 1:numel(e)
                bool_both = bldp_plot.flismember(e(i), e_both);
                bool_r = bldp_plot.flismember(e(i), e_r_filtered);
                bool_apx = bldp_plot.flismember(e(i), e_apx_filtered);
                if ~ (bool_both || bool_r || bool_apx)
                    e_filtered = [e_filtered e(i)];
                end
            end
            if numel(e_r_filtered) ~= numel(e_apx_filtered)
                error('We have numel(e_r_filtered) ~= numel(e_apx_filtered)');
            end
            
            if log_scale
                set(gca, 'YScale', 'log');
                ylim([10^(-4) 10^4])
            end

            if numel(e_r_filtered) > 0
                yline(t(e_r_filtered), 'LineStyle', '-', 'LineWidth', 1, 'Alpha', 0.25, 'Color', svd_colour); hold on;
            end
            if numel(e_apx_filtered) > 0
                yline(t(e_apx_filtered), 'LineStyle', ':', 'LineWidth', 2, 'Alpha', 0.25, 'Color', breg_colour); hold on;
            end
            
            %scatter(e, zeros(1, numel(e)), sz_original, 'black', 'filled', 'o', 'MarkerEdgeAlpha', 0.0, 'MarkerFaceAlpha', 0.8); hold on;
            s1 = scatter(e_apx_filtered, t(e_apx_filtered), sz_breg, 'filled', 'o', 'MarkerEdgeAlpha', .10, 'MarkerFaceAlpha', config.alpha, 'MarkerEdgeColor', breg_colour, 'MarkerFaceColor', breg_colour); hold on;
            s2 = scatter(e_r_filtered, t(e_r_filtered), sz_svd, 'filled', 'd', 'MarkerEdgeAlpha', 1.0, 'MarkerFaceAlpha', config.alpha, 'MarkerEdgeColor', '#00ecec', 'MarkerFaceColor', svd_colour); hold on;
            s3 = scatter(e_filtered, t(e_filtered), sz_original, 'black', 'filled', 'o', 'MarkerEdgeAlpha', 0.0, 'MarkerFaceAlpha', config.alpha); hold on;
            s4 = scatter(e_both, t(e_both), sz_both, 'filled', 'square', 'MarkerEdgeAlpha', 0.0, 'MarkerFaceAlpha', config.alpha, 'MarkerEdgeColor', config.both_colour, 'MarkerFaceColor', config.both_colour); hold on;
            
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
            ylabel('$\lambda-\log(1+\lambda)$', 'Interpreter', 'latex', 'FontSize', font_size)
            axis square;
            hold off;
            ldg = legend([s2 s1 s4], legends);
            set(ldg, 'Interpreter', 'latex', 'FontSize', font_size, 'Location', 'northoutside', 'Orientation', 'horizontal', 'NumColumns', 4);
            grid off;
            filename = 'bregman_curve.pdf';
            exportgraphics(gcf, fullfile(path, filename));
        end
      
        function plot_svd_curve(e, e_r, e_apx, path, ylog)
            font_size = 16;
            config = Plotting();
            svd_colour = config.svd.colour;
            breg_colour = config.breg.colour;

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
                if bldp_plot.flismember(e_r(i), e_apx)
                    if abs(e_r(i)) > 1e-10
                        e_both = [e_both e_r(i)];
                    end
                else
                    e_r_filtered = [e_r_filtered e_r(i)];
                end
            end
            for i = 1:r
                if ~ bldp_plot.flismember(e_apx(i), e_both)
                    e_apx_filtered = [e_apx_filtered e_apx(i)];
                end
            end

            for i = 1:numel(e)
                bool_both = bldp_plot.flismember(e(i), e_both);
                bool_r = bldp_plot.flismember(e(i), e_r_filtered);
                bool_apx = bldp_plot.flismember(e(i), e_apx_filtered);
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

            if ylog
                set(gca, 'YScale', 'log');
                ylim([10^(-4) 10^4])
            end


            if numel(e_r_filtered) ~= 0
                yline(t(shf + e_r_filtered), 'LineStyle', '-', 'LineWidth', 1, 'Alpha', 0.25, 'Color', svd_colour); hold on;
            end
            if numel(e_apx_filtered) ~= 0
                yline(t(shf + e_apx_filtered), 'LineStyle', ':', 'LineWidth', 2, 'Alpha', 0.25, 'Color', breg_colour); hold on;
            end
            
            xs = linspace(min(e), max(e)+1e-01, 100); hold on;
            plot(xs, arrayfun(t, xs), 'black');
            
            e_zeros = zeros(1, numel(e));
            %scatter(e, e_zeros, sz_original, 'black', 'filled', 'o', 'MarkerEdgeAlpha', 0.0, 'MarkerFaceAlpha', 0.8); hold on;
            s1 = scatter(shf + e_apx_filtered, t(shf + e_apx_filtered), sz_breg, 'filled', 'o', 'MarkerEdgeAlpha', 0.0, 'MarkerFaceAlpha', config.alpha, 'MarkerEdgeColor', breg_colour, 'MarkerFaceColor', breg_colour); hold on;
            s2 = scatter(shf + e_r_filtered, t(shf + e_r_filtered), sz_svd, 'filled', 'd', 'MarkerEdgeAlpha', 1.0, 'MarkerFaceAlpha', config.alpha, 'MarkerEdgeColor', '#00ecec', 'MarkerFaceColor', svd_colour); hold on;
            s3 = scatter(shf + e_filtered, t(shf + e_filtered), sz_original, 'black', 'filled', 'o', 'MarkerEdgeAlpha', 0.0, 'MarkerFaceAlpha', config.alpha);
            s4 = scatter(shf + e_both, t(shf + e_both), sz_both, 'filled', 'square', 'MarkerEdgeAlpha', 0.0, 'MarkerFaceAlpha', config.alpha, 'MarkerEdgeColor', config.both_colour, 'MarkerFaceColor', config.both_colour); hold on;
            
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
            ldg = legend([s2 s1 s4], {'$\left[ \! \left[ G\right] \! \right]_r$ (TSVD)', '$\langle\!\langle G\rangle\!\rangle_r$ (BT)', 'Overlap'});
            set(ldg, 'Interpreter', 'latex', 'FontSize', font_size, 'Location', 'northoutside', 'Orientation', 'horizontal', 'NumColumns', 4);%, 'Location', 'northoutside', 'Orientation', 'horizontal', 'NumColumns', 4);
            grid off;
            filename = 'svd_curve.pdf';
            exportgraphics(gcf, fullfile(path, filename));
        end
       
        function b = flismember(x, xs)
            b = 0;
            for i = 1:numel(xs)
               if abs(x - xs(i)) < 1e-12
                   b = 1;
               end
            end
        end
    end
end