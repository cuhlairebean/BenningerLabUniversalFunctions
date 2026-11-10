function keep = removeCellsGUI(CellTC, cellLabels, titleStr, displayData, yLabelStr)
% REMOVECELLSGUI  Plot all cell traces; hover to see the cell name, click a trace to remove it.
%   keep = removeCellsGUI(CellTC, cellLabels, titleStr, displayData, yLabelStr)
%   CellTC      : frames x cells matrix
%   cellLabels  : cell array of names (one per column)
%   displayData : optional frames x cells matrix shown INSTEAD of CellTC (e.g. the normalized dF/F)
%   yLabelStr   : optional y-axis label for displayData
%   keep        : logical row vector (1 x cells), false = removed.
%
%   - HOVER over the plot: the nearest trace turns red and a tooltip shows its cell label.
%   - CLICK: asks "Remove <cell>?" (Yes removes it).
%   - ZOOM: mouse scroll wheel zooms around the cursor; "Reset zoom" button (or key R) resets;
%           the toolbar zoom/pan tools also work (click-to-remove is paused while a tool is on).
%   - DONE: press the Done button, ENTER, or close the window.
yLabelUse = 'Intensity';
if nargin >= 4 && ~isempty(displayData) && isequal(size(displayData),size(CellTC))
    CellTC = displayData; yLabelUse = 'Normalized \DeltaF/F';
end
if nargin >= 5 && ~isempty(yLabelStr), yLabelUse = yLabelStr; end
nC = size(CellTC,2);
nT = size(CellTC,1);
if nargin < 2 || isempty(cellLabels)
    cellLabels = arrayfun(@(j) sprintf('cell%d',j),1:nC,'UniformOutput',false);
end
if nargin < 3, titleStr = ''; end
keep = true(1,nC);
if nC == 0, return; end
x = (1:nT)';

fig = figure('Name','Remove outlier cells','NumberTitle','off');
ax = axes('Parent',fig,'Position',[0.08 0.14 0.89 0.76]); hold(ax,'on');
h = gobjects(1,nC);
for c = 1:nC
    h(c) = plot(ax,x,CellTC(:,c),'Color',[0.3 0.3 0.3 0.6],'HitTest','off');
end
xlabel(ax,'Frame index'); ylabel(ax,yLabelUse);
axis(ax,'tight');
xl0 = xlim(ax); yl0 = ylim(ax);
tip = text(ax,0,0,'','Visible','off','BackgroundColor',[1 1 0.85],'EdgeColor',[0.4 0.4 0.4], ...
    'Margin',3,'FontWeight','bold','Interpreter','none','VerticalAlignment','bottom','HitTest','off');
uicontrol(fig,'Style','pushbutton','String','Done','Units','normalized','Position',[0.88 0.01 0.1 0.06], ...
    'FontWeight','bold','Callback',@(~,~) uiresume(fig));
uicontrol(fig,'Style','pushbutton','String','Reset zoom','Units','normalized','Position',[0.76 0.01 0.11 0.06], ...
    'Callback',@(~,~) resetView());

hover = 0;
refreshTitle();
set(fig,'WindowButtonMotionFcn',@onMove,'WindowButtonDownFcn',@onClick, ...
    'WindowScrollWheelFcn',@onScroll,'KeyPressFcn',@onKey,'CloseRequestFcn',@(~,~) uiresume(fig));
uiwait(fig);
if isvalid(fig), delete(fig); end

    function refreshTitle()
        title(ax,sprintf('%s | %d of %d cells kept | hover = name, click = remove, scroll = zoom, ENTER/Done = finish', ...
            titleStr,nnz(keep),nC),'Interpreter','none');
    end

    function [c,ok] = nearestCell()
        ok = false; c = 0;
        cp = get(ax,'CurrentPoint'); px = cp(1,1); py = cp(1,2);
        xl = xlim(ax); yl = ylim(ax);
        if px < xl(1) || px > xl(2) || py < yl(1) || py > yl(2), return; end
        fr = min(max(round(px),1),nT);
        d = abs(CellTC(fr,:) - py)/max(diff(yl),eps);
        d(~keep) = inf;
        [dmin,c] = min(d);
        ok = isfinite(dmin);
    end

    function setHover(c)
        if hover > 0 && hover <= nC && isvalid(h(hover)) && keep(hover)
            set(h(hover),'Color',[0.3 0.3 0.3 0.6],'LineWidth',0.5);
        end
        hover = c;
        if c > 0
            set(h(c),'Color',[1 0 0],'LineWidth',2);
            uistack(h(c),'top');
        end
    end

    function onMove(~,~)
        if ~isvalid(fig), return; end
        [c,ok] = nearestCell();
        if ~ok
            setHover(0); set(tip,'Visible','off'); return;
        end
        if c ~= hover, setHover(c); end
        cp = get(ax,'CurrentPoint');
        set(tip,'String',cellLabels{c},'Position',[cp(1,1) cp(1,2)],'Visible','on');
        uistack(tip,'top');
    end

    function onClick(~,~)
        if ~strcmp(get(fig,'SelectionType'),'normal'), return; end
        [c,ok] = nearestCell();
        if ~ok, return; end
        setHover(c);
        ans_ = questdlg(sprintf('Remove %s?',cellLabels{c}),'Remove cell','Yes','No','No');
        if strcmp(ans_,'Yes')
            keep(c) = false;
            set(h(c),'Visible','off');
            hover = 0; set(tip,'Visible','off');
        else
            setHover(0);
        end
        refreshTitle();
    end

    function onScroll(~,evt)
        cp = get(ax,'CurrentPoint'); px = cp(1,1); py = cp(1,2);
        xl = xlim(ax); yl = ylim(ax);
        if px < xl(1) || px > xl(2) || py < yl(1) || py > yl(2), return; end
        f = 1.25^evt.VerticalScrollCount;          % scroll up (negative) = zoom in
        nx = px + (xl - px)*f; ny = py + (yl - py)*f;
        nx(1) = max(nx(1),xl0(1)); nx(2) = min(nx(2),xl0(2));
        ny(1) = max(ny(1),yl0(1)-0.5*diff(yl0)); ny(2) = min(ny(2),yl0(2)+0.5*diff(yl0));
        if nx(2) - nx(1) > 2 && ny(2) > ny(1), xlim(ax,nx); ylim(ax,ny); end
    end

    function onKey(~,evt)
        switch lower(evt.Key)
            case {'return','enter'}, uiresume(fig);
            case 'r', resetView();
        end
    end

    function resetView()
        xlim(ax,xl0); ylim(ax,yl0);
    end
end
