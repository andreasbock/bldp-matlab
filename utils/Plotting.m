function config = Plotting()

    alpha = 0.3;
    size_step = 50;
    line_width = 1.8;

    config.font_size = 16;
    config.font_size_legend = 12;
    config.font_size_axes = 20;
    config.font_size_title = 20;
    config.font_size_ticks = 14;

    config.both_colour = '#379c37';
    config.alpha = alpha;

    %% No preconditioner
    config.nopc.line_width = line_width;
    config.nopc.line_style = ':';
    config.nopc.label = '$I$';
    config.nopc.marker = '*';
    config.nopc.colour = 'red';
    config.nopc.alpha = alpha;
    config.nopc.size_step = size_step;

    %% Incomplete Cholesky
    config.ichol.line_width = line_width;
    config.ichol.line_style = '-';
    config.ichol.label = '$Q$';
    config.ichol.marker = 'square';
    config.ichol.colour = 'blue';
    config.ichol.alpha = alpha;
    config.ichol.size_step = size_step;

    %% SVD
    config.svd.line_width = line_width;
    config.svd.line_style = '--';
    config.svd.label = 'TSVD';
    config.svd.marker = 'd';
    config.svd.colour = 'cyan';
    config.svd.alpha = alpha;
    config.svd.size_step = size_step;

    %% Nyström
    config.nystrom.line_width = line_width;
    config.nystrom.line_style = '-.';
    config.nystrom.label = 'Nyström';
    config.nystrom.marker = 'd';
    config.nystrom.colour = 'cyan';
    config.nystrom.alpha = alpha;
    config.nystrom.size_step = size_step;

    %% Bregman
    config.breg.line_width = line_width;
    config.breg.line_style = ':';
    config.breg.label = 'Bregman';
    config.breg.marker = 'square';
    config.breg.colour = 'magenta';
    config.breg.alpha = alpha;
    config.breg.size_step = size_step;

    %% Reverse Bregman
    config.rbreg.line_width = line_width;
    config.rbreg.line_style = '-';
    config.rbreg.label = 'Rev. Bregman';
    config.rbreg.marker = 'pentagram';
    config.rbreg.colour = '#8e61c2';
    config.rbreg.alpha = alpha;
    config.rbreg.size_step = size_step;

    font_size_step = 25;
    step = 10;
end

