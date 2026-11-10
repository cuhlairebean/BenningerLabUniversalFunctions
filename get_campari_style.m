function style = get_campari_style(n_conds)
% GET_CAMPARI_STYLE Returns a struct containing styling parameters for the CaMPARI theme.
    if nargin < 1, n_conds = 3; end
    style.colors = containers.Map(...
        {'CaMPARI green', 'CaMPARI Red', 'green', 'magenta', 'Green', 'Yellow', 'Red'}, ...
        {[0.15, 0.75, 0.15], [0.75, 0.10, 0.75], [0.15, 0.75, 0.15], [0.75, 0.10, 0.75], ...
         [0.15, 0.75, 0.15], [0.85, 0.75, 0.10], [0.75, 0.10, 0.75]});
        
    style.hatch_colors = containers.Map(...
        {'CaMPARI green', 'CaMPARI Red', 'green', 'magenta', 'Green', 'Yellow', 'Red'}, ...
        {[0.05, 0.50, 0.05], [0.50, 0.05, 0.50], [0.05, 0.50, 0.05], [0.50, 0.05, 0.50], ...
         [0.05, 0.50, 0.05], [0.65, 0.55, 0.05], [0.50, 0.05, 0.50]});
        
    raw_patterns = {'none', 'horizontal', 'diag_left', 'diag_right'};
    style.cond_hatch = repmat(raw_patterns, 1, ceil(n_conds / 4));
    style.cond_hatch = style.cond_hatch(1:n_conds);
    
    style.font_name        = 'Avenir';
    style.title_size       = 16;
    style.axis_size        = 15;
    style.box_width        = 0.65;
    style.num_hatch_lines  = 5;
    style.hatch_line_width = 1.0;
    style.fig_width        = 750;
    style.fig_height       = 650;
    style.face_alpha       = 0.45;
end
