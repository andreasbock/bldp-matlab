function lh = plot_resvec(resvec_long, config)
    last_indx = find(resvec_long > 0, 1, 'last');
    resvec = resvec_long(1:last_indx);
    lh_ = semilogy(0:length(resvec)-1, resvec, 'DisplayName', config.label, 'LineStyle', config.line_style, 'LineWidth', config.line_width, 'Color', config.colour, 'MarkerEdgeColor', config.colour, 'HandleVisibility', 'off'); hold on;
    lh_.Color = [lh_.Color config.alpha];
    hold on;
    lh = scatter(0:length(resvec)-1, resvec, config.size_step, config.marker, 'DisplayName', config.label, 'MarkerFaceColor', config.colour, 'MarkerEdgeColor', config.colour, 'MarkerFaceAlpha', config.alpha, 'MarkerEdgeAlpha', config.alpha);
    hold on;
end