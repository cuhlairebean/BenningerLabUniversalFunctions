function keep = removeCellsGUI(CellTC, cellLabels, titleStr)
% REMOVECELLSGUI  Plot all cell traces; click a trace to be asked "Remove cell x?".
%   keep = removeCellsGUI(CellTC, cellLabels, titleStr)
%   CellTC     : frames x cells matrix
%   cellLabels : cell array of names (one per column)
%   keep       : logical row vector (1 x cells), false = removed.
% Click near a trace -> the nearest trace is highlighted and you are asked
% whether to remove it. Press ENTER (or right-click) in the plot when done.
nC = size(CellTC,2);
if nargin < 2 || isempty(cellLabels)
    cellLabels = arrayfun(@(j) sprintf('cell%d',j),1:nC,'UniformOutput',false);
end
if nargin < 3, titleStr = ''; end
keep = true(1,nC);
x = (1:size(CellTC,1))';

fig = figure('Name','Remove outlier cells','NumberTitle','off');
ax = axes('Parent',fig); hold(ax,'on');
h = gobjects(1,nC);
for c = 1:nC
    h(c) = plot(ax,x,CellTC(:,c),'Color',[0.3 0.3 0.3 0.6]);
end
xlabel(ax,'Frame index'); ylabel(ax,'Intensity');
done = false;
while ~done && isvalid(fig)
    title(ax,sprintf('%s | %d of %d cells kept | click a trace to remove it, ENTER = done', ...
        titleStr,nnz(keep),nC),'Interpreter','none');
    figure(fig);
    [px,py,btn] = ginput(1);
    if isempty(px) || btn ~= 1, break; end
    % nearest kept trace to the click (distance at the clicked frame, in y)
    fr = min(max(round(px),1),size(CellTC,1));
    d = abs(CellTC(fr,:) - py);
    d(~keep) = inf;
    [dmin,c] = min(d);
    if ~isfinite(dmin), continue; end
    set(h(c),'Color',[1 0 0],'LineWidth',2);
    uistack(h(c),'top');
    drawnow;
    ans_ = questdlg(sprintf('Remove %s?',cellLabels{c}),'Remove cell','Yes','No','No');
    if strcmp(ans_,'Yes')
        keep(c) = false;
        set(h(c),'Visible','off');
    else
        set(h(c),'Color',[0.3 0.3 0.3 0.6],'LineWidth',0.5);
    end
end
if isvalid(fig), close(fig); end
end
