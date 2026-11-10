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

ALPHA = 0.3;
nC = size(N,2);
campari = nargin >= 6 && ~isempty(ratio) && ~isnan(threshold) && numel(ratio) == nC;
fig = figure('Name',['Normalized traces ' ttl]);
ax = axes('Parent',fig); hold(ax,'on');
cls = repmat({'all'},1,nC);
means.all = mean(N,2,'omitnan'); means.red = nan(size(means.all)); means.green = nan(size(means.all));

if ~campari
    for c = 1:nC
        plot(ax,t,N(:,c),'Color',[0.5 0.5 0.5 ALPHA],'LineWidth',0.5);
    end
    hAll = plot(ax,t,means.all,'k','LineWidth',2.5);
    legend(ax,hAll,sprintf('Mean, all cells (n=%d)',nC),'Location','best','Box','off');
else
    ratio = ratio(:)';
    gCol = [0 0.75 0]; mCol = [1 0 1];
    dev = max(abs(ratio(~isnan(ratio)) - threshold));
    if isempty(dev) || dev == 0, dev = 1; end
    isRed = ratio >= threshold; isGreen = ratio < threshold;      % NaN -> both false
    cls(isRed) = {'red'}; cls(isGreen) = {'green'}; cls(isnan(ratio)) = {'unclassified'};
    for c = 1:nC
        if isnan(ratio(c)), col = [0.5 0.5 0.5];
        else
            f = min(max(0.5 + 0.5*(ratio(c)-threshold)/dev, 0), 1);   % 0.5 = at threshold
            col = (1-f)*gCol + f*mCol;
        end
        plot(ax,t,N(:,c),'Color',[col ALPHA],'LineWidth',0.5);
    end
    if any(isRed),   means.red   = mean(N(:,isRed),2,'omitnan');   end
    if any(isGreen), means.green = mean(N(:,isGreen),2,'omitnan'); end
    hAll = plot(ax,t,means.all,'k','LineWidth',2.5);
    hR = []; hG = []; lbl = {sprintf('Mean, all cells (n=%d)',nC)}; hh = hAll;
    if any(isRed),   hR = plot(ax,t,means.red,'Color',mCol,'LineWidth',2.5);  hh(end+1) = hR; lbl{end+1} = sprintf('Mean, red cells (n=%d)',nnz(isRed)); end
    if any(isGreen), hG = plot(ax,t,means.green,'Color',gCol,'LineWidth',2.5); hh(end+1) = hG; lbl{end+1} = sprintf('Mean, green cells (n=%d)',nnz(isGreen)); end
    legend(ax,hh,lbl,'Location','best','Box','off');
    % colour scale = the red:green ratio, threshold in the middle
    nCol = 64; f = linspace(0,1,nCol)';
    colormap(ax,(1-f)*gCol + f*mCol);
    caxis(ax,[threshold-dev threshold+dev]);
    cb = colorbar(ax); cb.Label.String = sprintf('Red:Green ratio (threshold = %.3g)',threshold);
end
xlabel(ax,xlab); ylabel(ax,'\DeltaF/F');
title(ax,ttl,'Interpreter','none');
box(ax,'off');
end
