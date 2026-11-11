function [fig, cls, means] = plotNormalizedTraces(t, N, xlab, ttl, ratio, threshold)
%PLOTNORMALIZEDTRACES  Normalized (dF/F) traces of one islet: every cell semi-transparent, bold opaque means.
%
%   [fig, cls, means] = plotNormalizedTraces(t, N, xlab, ttl)                       plain figure
%   [fig, cls, means] = plotNormalizedTraces(t, N, xlab, ttl, ratio, threshold)     CaMPARI figure
%
%   t          x values (time in s or frame index), N = frames x cells normalized matrix
%   ratio      1 x cells red:green ratio of each plotted cell (NaN = unknown), threshold = red/green cut-off
%   PLAIN:   all cells grey (alpha 0.3) + bold BLACK mean of all cells.
%   CAMPARI: each cell coloured by its ratio (green = low -> magenta = high, centred on the threshold,
%            alpha 0.3) + bold opaque means: BLACK = all cells, MAGENTA = red cells (ratio >= threshold),
%            GREEN = green cells (ratio < threshold). Cells with NaN ratio are grey and in no colour mean.
%   cls        cellstr per cell: 'red' / 'green' / 'unclassified' (all 'all' for the plain figure)
%   means      struct with .all .red .green (frames x 1, NaN if that group is empty)

ALPHA = 0.1;
nC = size(N,2);
campari = nargin >= 6 && ~isempty(ratio) && ~isnan(threshold) && numel(ratio) == nC;
% CaMPARI theme (colors/font/sizes) is used ONLY for CaMPARI figures; the main script passes ratio/threshold
% only when the user is Evelyn, so everyone else gets the plain default figure.
if campari, sty = get_campari_style(); end
fig = figure('Name',['Normalized traces ' ttl]);
ax = axes('Parent',fig); hold(ax,'on');
cls = repmat({'all'},1,nC);
means.all = mean(N,2,'omitnan'); means.red = nan(size(means.all)); means.green = nan(size(means.all));

if ~campari
    for c = 1:nC
        plot(ax,t,N(:,c),'Color',[0.5 0.5 0.5 ALPHA],'LineWidth',0.5);
    end
    hAll = plot(ax,t,means.all,'k','LineWidth',2.5);
    legend(ax,hAll,sprintf('Mean, all cells (n=%d)',nC),'Location','southoutside','Orientation','horizontal','Box','off','FontSize',9,'TextColor','k');
else
    ratio = ratio(:)';
    gCol = sty.colors('green'); mCol = sty.colors('magenta'); yCol = sty.colors('Yellow');
    % gradient by colour stops (position 0-1 along the colorbar, 0.5 = threshold), linear between stops:
    %   0% green #339933 | 34% #ADDE43 | 50% #D9BF1A (threshold) | 62% #E0B9ED | 100% magenta #D784FF
    stopPos = [0 0.34 0.50 0.62 1];
    stopCol = [gCol; [173 222 67]/255; [217 191 26]/255; [224 185 237]/255; mCol];
    grad = @(f) interp1(stopPos,stopCol,min(max(f,0),1));
    % Yellow (#D9BF1A) sits at ratio = threshold ALWAYS (0 on the "ratio - threshold" axis). Each side has its own
    % scale: below the threshold the green end is the lowest ratio, above it the magenta end is the highest ratio,
    % and the intermediate stops (34% / 62%) are placed relative to that, so yellow never moves.
    rv = ratio(~isnan(ratio));
    dmin = max(0,threshold - min(rv)); dmax = max(0,max(rv) - threshold);
    if isempty(rv) || (dmin == 0 && dmax == 0), dmin = 1; dmax = 1; end
    if dmin == 0, dmin = dmax; end
    if dmax == 0, dmax = dmin; end
    fmap = @(r) min(max(0.5 + 0.5*(r-threshold)./(dmax*(r>=threshold) + dmin*(r<threshold)),0),1);
    isRed = ratio >= threshold; isGreen = ratio < threshold;      % NaN -> both false
    cls(isRed) = {'red'}; cls(isGreen) = {'green'}; cls(isnan(ratio)) = {'unclassified'};
    for c = 1:nC
        if isnan(ratio(c)), col = [0.5 0.5 0.5];
        else
            col = grad(fmap(ratio(c)));
        end
        plot(ax,t,N(:,c),'Color',[col ALPHA],'LineWidth',0.5);
    end
    if any(isRed),   means.red   = mean(N(:,isRed),2,'omitnan');   end
    if any(isGreen), means.green = mean(N(:,isGreen),2,'omitnan'); end
    hAll = plot(ax,t,means.all,'k','LineWidth',2.5);
    hR = []; hG = []; lbl = {sprintf('Mean, all cells (n=%d)',nC)}; hh = hAll;
    if any(isRed),   hR = plot(ax,t,means.red,'Color',mCol,'LineWidth',2.5);  hh(end+1) = hR; lbl{end+1} = sprintf('Mean, red cells (n=%d)',nnz(isRed)); end
    if any(isGreen), hG = plot(ax,t,means.green,'Color',gCol,'LineWidth',2.5); hh(end+1) = hG; lbl{end+1} = sprintf('Mean, green cells (n=%d)',nnz(isGreen)); end
    legend(ax,hh,lbl,'Location','southoutside','Orientation','horizontal','Box','off','FontSize',9,'TextColor','k');
    % colour scale = the red:green ratio; the colormap follows the same mapping, so yellow is at the threshold
    nCol = 128; vals = linspace(threshold-dmin,threshold+dmax,nCol)';
    cm = zeros(nCol,3);
    for q = 1:nCol, cm(q,:) = grad(fmap(vals(q))); end
    colormap(ax,cm);
    caxis(ax,[threshold-dmin threshold+dmax]);
    cb = colorbar(ax); cb.Label.String = sprintf('Red:Green ratio (threshold = %.2f)',threshold);
    % ticks with 2 decimals; ticks too close to the threshold are dropped so its label stays readable
    lo = threshold-dmin; hi = threshold+dmax;
    tk = cb.Ticks; tk = tk(abs(tk-threshold) > 0.08*(hi-lo));
    cb.Ticks = sort([tk, threshold]);
    cb.TickLabels = arrayfun(@(v) sprintf('%.2f',v),cb.Ticks,'UniformOutput',false);
    % the threshold gets a grey, translucent, bold mark just across the colorbar
    thrLine = annotation(fig,'line',[0 0],[0 0],'LineWidth',3);
    try, thrLine.Color = [0.35 0.35 0.35 0.6]; catch, thrLine.Color = [0.6 0.6 0.6]; end   % alpha needs a recent MATLAB
    set(fig,'SizeChangedFcn',@(~,~) placeThrMark(cb,thrLine,threshold,lo,hi,ax));
    drawnow; placeThrMark(cb,thrLine,threshold,lo,hi,ax);
end
xlabel(ax,xlab); ylabel(ax,'\DeltaF/F');
title(ax,ttl,'Interpreter','none');
box(ax,'off');
if campari   % CaMPARI theme font (falls back to the default font if not installed)
    fn = sty.font_name; if ~any(strcmp(listfonts,fn)), fn = get(groot,'defaultAxesFontName'); end
    set(ax,'FontName',fn,'FontSize',sty.axis_size);
    ax.Title.FontSize = sty.title_size; ax.Title.FontName = fn;
    ax.XLabel.FontSize = sty.axis_size; ax.YLabel.FontSize = sty.axis_size;
end
applyLightTheme(fig);   % white background, black fonts/axes
if campari, drawnow; placeThrMark(cb,thrLine,threshold,lo,hi,ax); end   % layout may have shifted with the fonts
end

function placeThrMark(cb,ln,thr,lo,hi,ax)
% 1) colorbar height = 25% of the plot height, vertically centred on the plot;
% 2) short grey line across the colorbar at the threshold (figure-normalized units; follows on resize)
if ~isvalid(cb) || ~isvalid(ln), return; end
oldU = cb.Units; cb.Units = 'normalized';
oldA = ax.Units; ax.Units = 'normalized'; ap = ax.Position; ax.Units = oldA;
p = cb.Position;
p(4) = 0.25*ap(4); p(2) = ap(2) + 0.5*(ap(4)-p(4));
cb.Position = p; cb.Units = oldU;
y = p(2) + p(4)*(thr-lo)/(hi-lo);
set(ln,'X',[p(1)-0.15*p(3), p(1)+1.15*p(3)],'Y',[y y]);
end
