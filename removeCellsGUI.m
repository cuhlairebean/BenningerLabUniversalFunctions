function keep = removeCellsGUI(CellTC, cellLabels, titleStr, displayData, yLabelStr, maskImg, colIDs)
% REMOVECELLSGUI  Plot all cell traces; hover to see the cell name, click a trace to remove it.
%   keep = removeCellsGUI(CellTC, cellLabels, titleStr, displayData, yLabelStr, maskImg, colIDs)
%   CellTC      : frames x cells matrix
%   cellLabels  : cell array of names (one per column)
%   displayData : optional frames x cells matrix shown INSTEAD of CellTC (e.g. the normalized dF/F)
%   yLabelStr   : optional y-axis label for displayData
%   maskImg     : optional label image of the cells (pixel value = cell ID, 0 = background)
%   colIDs      : cell ID (value in maskImg) of each column of CellTC
%   keep        : logical row vector (1 x cells), false = removed.
%
%   With maskImg/colIDs a second panel shows the CELL MASKS. A cell can then be removed EITHER way:
%     - click its TRACE in the time course, or
%     - click the CELL in the mask panel.
%   Hovering either one highlights the cell (red) in both panels.
%
%   - HOVER: the nearest trace / the cell under the cursor turns red and a tooltip shows its label.
%   - CLICK: asks "Remove <cell>?" (Yes removes it).
%   - ZOOM: mouse scroll wheel zooms around the cursor (either panel); "Reset zoom" button (or key R) resets;
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
hasMask = nargin >= 7 && ~isempty(maskImg) && numel(colIDs) == nC;
keep = true(1,nC);
if nC == 0, return; end
x = (1:nT)';

fig = figure('Name','Remove outlier cells','NumberTitle','off');
if hasMask
    set(fig,'Units','normalized','Position',[0.05 0.1 0.9 0.75]);
    ax = axes('Parent',fig,'Position',[0.05 0.14 0.50 0.76]);
else
    ax = axes('Parent',fig,'Position',[0.08 0.14 0.89 0.76]);
end
hold(ax,'on');
h = gobjects(1,nC);
for ci = 1:nC
    h(ci) = plot(ax,x,CellTC(:,ci),'Color',[0.3 0.3 0.3 0.6],'HitTest','off');
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

% --- optional mask panel: lut maps a pixel value (cell ID) to the column of that cell ---
axM = []; hMask = []; colImg = []; tipM = []; xlM0 = []; ylM0 = [];
if hasMask
    maskImg = round(double(maskImg));
    maskImg(maskImg < 0 | isnan(maskImg)) = 0;
    lut = zeros(1,max([maskImg(:); colIDs(:)])+1);
    lut(round(colIDs(:))'+1) = 1:nC;
    colImg = reshape(lut(maskImg+1),size(maskImg));        % 0 = background, else column index
    axM = axes('Parent',fig,'Position',[0.58 0.14 0.40 0.76]);
    hMask = image(axM,renderMask(0)); axis(axM,'image'); axis(axM,'off');
    set(hMask,'HitTest','off'); set(axM,'HitTest','off');
    title(axM,'Cell masks: click a cell to remove it','Interpreter','none');
    hold(axM,'on');
    tipM = text(axM,0,0,'','Visible','off','BackgroundColor',[1 1 0.85],'EdgeColor',[0.4 0.4 0.4], ...
        'Margin',3,'FontWeight','bold','Interpreter','none','VerticalAlignment','bottom','HitTest','off');
    xlM0 = xlim(axM); ylM0 = ylim(axM);
end

hover = 0; busy = false;
refreshTitle();
set(fig,'WindowButtonMotionFcn',@onMove,'WindowButtonDownFcn',@onClick, ...
    'WindowScrollWheelFcn',@onScroll,'KeyPressFcn',@onKey,'CloseRequestFcn',@(~,~) uiresume(fig));
uiwait(fig);
if isvalid(fig), delete(fig); end

    function im = renderMask(hv)
        % background dark, kept cells light grey, hovered cell red (removed cells disappear)
        im = zeros([size(colImg) 3]);
        kp = find(keep);
        k1 = ismember(colImg,kp) & colImg > 0;
        for ch = 1:3, tmp = im(:,:,ch); tmp(k1) = 0.65; im(:,:,ch) = tmp; end
        if hv > 0
            k2 = colImg == hv;
            im(:,:,1) = im(:,:,1) + 0*k2; tmp = im(:,:,1); tmp(k2) = 1; im(:,:,1) = tmp;
            tmp = im(:,:,2); tmp(k2) = 0; im(:,:,2) = tmp;
            tmp = im(:,:,3); tmp(k2) = 0; im(:,:,3) = tmp;
        end
    end

    function refreshTitle()
        if hasMask
            how = 'click a trace OR a cell in the mask to remove';
        else
            how = 'click a trace to remove';
        end
        title(ax,sprintf('%s | %d of %d cells kept | hover = name | %s | scroll = zoom | ENTER/Done = finish', ...
            titleStr,nnz(keep),nC,how),'Interpreter','none');
    end

    function r = inAxes(a)
        r = false;
        if isempty(a) || ~isvalid(a), return; end
        cp = get(a,'CurrentPoint'); px = cp(1,1); py = cp(1,2);
        xl = xlim(a); yl = ylim(a);
        r = px >= xl(1) && px <= xl(2) && py >= yl(1) && py <= yl(2);
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

    function [c,ok] = maskCell()
        ok = false; c = 0;
        if ~hasMask || ~inAxes(axM), return; end
        cp = get(axM,'CurrentPoint');
        col = round(cp(1,1)); row = round(cp(1,2));
        if row < 1 || col < 1 || row > size(colImg,1) || col > size(colImg,2), return; end
        c = colImg(row,col);
        ok = c > 0 && keep(c);
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
        if hasMask, set(hMask,'CData',renderMask(c)); end
    end

    function onMove(~,~)
        if ~isvalid(fig) || busy, return; end
        set(tip,'Visible','off');
        if hasMask, set(tipM,'Visible','off'); end
        [c,ok] = maskCell();
        if ok
            if c ~= hover, setHover(c); end
            cp = get(axM,'CurrentPoint');
            set(tipM,'String',cellLabels{c},'Position',[cp(1,1) cp(1,2)],'Visible','on');
            return;
        end
        [c,ok] = nearestCell();
        if ~ok
            if hover ~= 0, setHover(0); end
            return;
        end
        if c ~= hover, setHover(c); end
        cp = get(ax,'CurrentPoint');
        set(tip,'String',cellLabels{c},'Position',[cp(1,1) cp(1,2)],'Visible','on');
        uistack(tip,'top');
    end

    function onClick(~,~)
        if busy, return; end                         % ignore clicks/hover while the dialog is open
        if ~strcmp(get(fig,'SelectionType'),'normal'), return; end
        [sel,okSel] = maskCell();                  % clicked a cell in the mask panel?
        if ~okSel, [sel,okSel] = nearestCell(); end   % else the nearest trace
        if ~okSel || sel < 1, return; end
        busy = true;                               % freeze: the cell is fixed to 'sel' from here on
        setHover(sel);
        try
            ans_ = questdlg(sprintf('Remove %s?',cellLabels{sel}),'Remove cell','Yes','No','No');
        catch
            ans_ = 'No';
        end
        if strcmp(ans_,'Yes')
            keep(sel) = false;
            set(h(sel),'Visible','off');
            hover = 0; set(tip,'Visible','off');
            if hasMask, set(hMask,'CData',renderMask(0)); set(tipM,'Visible','off'); end
        else
            setHover(0);
        end
        busy = false;
        refreshTitle();
    end

    function onScroll(~,evt)
        f = 1.25^evt.VerticalScrollCount;          % scroll up (negative) = zoom in
        if inAxes(ax)
            cp = get(ax,'CurrentPoint'); px = cp(1,1); py = cp(1,2);
            xl = xlim(ax); yl = ylim(ax);
            nx = px + (xl - px)*f; ny = py + (yl - py)*f;
            nx(1) = max(nx(1),xl0(1)); nx(2) = min(nx(2),xl0(2));
            ny(1) = max(ny(1),yl0(1)-0.5*diff(yl0)); ny(2) = min(ny(2),yl0(2)+0.5*diff(yl0));
            if nx(2) - nx(1) > 2 && ny(2) > ny(1), xlim(ax,nx); ylim(ax,ny); end
        elseif hasMask && inAxes(axM)
            cp = get(axM,'CurrentPoint'); px = cp(1,1); py = cp(1,2);
            xl = xlim(axM); yl = ylim(axM);
            nx = px + (xl - px)*f; ny = py + (yl - py)*f;
            nx(1) = max(nx(1),xlM0(1)); nx(2) = min(nx(2),xlM0(2));
            ny(1) = max(ny(1),ylM0(1)); ny(2) = min(ny(2),ylM0(2));
            if nx(2) - nx(1) > 4 && ny(2) - ny(1) > 4, xlim(axM,nx); ylim(axM,ny); end
        end
    end

    function onKey(~,evt)
        switch lower(evt.Key)
            case {'return','enter'}, uiresume(fig);
            case 'r', resetView();
        end
    end

    function resetView()
        xlim(ax,xl0); ylim(ax,yl0);
        if hasMask, xlim(axM,xlM0); ylim(axM,ylM0); end
    end
end
