function applyLightTheme(fig)
%APPLYLIGHTTHEME  Force the light theme on a figure: white background, black fonts, axes, ticks and lines.
%   applyLightTheme(fig)       fix an existing figure (axes, legend, colorbar, text)
%   applyLightTheme()          set the session defaults so every NEW figure is light
if nargin < 1
    set(groot,'defaultFigureColor','w','defaultAxesColor','w','defaultAxesXColor','k', ...
        'defaultAxesYColor','k','defaultAxesZColor','k','defaultTextColor','k', ...
        'defaultAxesGridColor','k','defaultAxesMinorGridColor','k');
    try, set(groot,'defaultLegendTextColor','k','defaultLegendColor','w','defaultLegendEdgeColor','k'); catch, end
    try, s = settings; s.matlab.appearance.MATLABTheme.TemporaryValue = 'light'; catch, end
    return
end
try, theme(fig,'light'); catch, end            % R2025a+: switch the whole figure to the light theme
set(fig,'Color','w');
ax = findall(fig,'Type','axes');
for a = ax(:)'
    set(a,'Color','w','XColor','k','YColor','k','ZColor','k','GridColor','k','MinorGridColor','k');
    try, a.Title.Color = 'k'; a.XLabel.Color = 'k'; a.YLabel.Color = 'k'; catch, end
end
lg = findall(fig,'Type','legend');
for l = lg(:)', try, set(l,'TextColor','k','Color','w','EdgeColor','k'); catch, end, end
cb = findall(fig,'Type','colorbar');
for c = cb(:)', try, set(c,'Color','k'); c.Label.Color = 'k'; catch, end, end
tx = findall(fig,'Type','text');
for t = tx(:)', try, if ~strcmp(get(t,'Visible'),'off') && isequal(get(t,'BackgroundColor'),'none'), set(t,'Color','k'); end, catch, end, end
end
