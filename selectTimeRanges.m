function [periods, tbl] = selectTimeRanges(tc, names, figTitle, fig, initial, pickMask)
%SELECTTIMERANGES  Click start/end of named periods on a timecourse plot.
%
%   [periods, tbl] = selectTimeRanges(tc, names, figTitle, fig, initial, pickMask)
%
%   initial   (optional) Nx2 [start end] of previously saved selections (NaN = not used).
%             If given, they are drawn and you go straight to the review menu, where you can
%             accept them or re-select any period.
%
%   pickMask  (optional) logical 1xN, used together with 'initial': only the periods marked true
%             are asked for in the first pass; the rest are taken from 'initial' as they are
%             (e.g. low/high/drug glucose already chosen for the first islet of the same file).
%
%   Every period can be entered by CLICKING START/END on the plot or by TYPING the start and
%   end frame indices (you choose the method; it is asked again when you re-select a period).
%
%   tc        vector (mean timecourse, e.g. one islet)
%   names     cellstr, e.g. {'Low glucose','High glucose','Drug','Baseline','First phase','Second phase'}
%   figTitle  title text
%   fig       (optional) existing figure holding the plot of tc; one is created if omitted
%
%   For every period you click START then END on the plot (overlap between periods is
%   fine). Afterwards a menu lets you ACCEPT, re-select any single period, or redo all.
%   Press Enter in the plot without clicking to mark a period as "not used".
%
%   periods   struct array: .name .start .end .indices   (NaN start/end = not used)
%   tbl       matrix, one column per period holding its frame indices
%             (shorter columns padded with NaN)

tc = tc(:);
n = numel(tc);
N = numel(names);
if nargin < 3, figTitle = ''; end
if nargin < 4 || isempty(fig) || ~isvalid(fig)
    fig = figure('Name',figTitle);
    plot(tc);
end
figure(fig); ax = gca; hold(ax,'on');
colors = lines(N);
st = nan(1,N); en = nan(1,N);
hP = gobjects(1,N); hT = gobjects(1,N);

if nargin >= 5 && ~isempty(initial) && size(initial,1) == N
    for i = 1:N
        if ~any(isnan(initial(i,:)))
            st(i) = max(1,round(min(initial(i,:)))); en(i) = min(n,round(max(initial(i,:))));
        end
        drawPatch(i);
    end
    if nargin >= 6 && ~isempty(pickMask)
        for i = find(pickMask(:)'), pickRange(i); end
    end
else
    for i = 1:N, pickRange(i); end
end

while true
    items = cell(1,N+2);
    items{1} = '>>> ACCEPT these selections <<<';
    for i = 1:N
        if isnan(st(i)), r = 'not used'; else, r = sprintf('%d - %d',st(i),en(i)); end
        items{i+1} = sprintf('Re-select: %s  [%s]',names{i},r);
    end
    items{N+2} = 'Redo ALL periods';
    [sel,ok] = listdlg('PromptString',[figTitle ' - check ranges:'],'ListString',items, ...
        'SelectionMode','single','ListSize',[380 200],'Name','Period selection');
    if ~ok || sel == 1, break; end          % closing the dialog also accepts
    if sel == N+2
        for i = 1:N, pickRange(i); end
    else
        pickRange(sel-1);
    end
end
title(ax,figTitle,'Interpreter','none');

periods = struct('name',names(:),'start',num2cell(st(:)),'end',num2cell(en(:)),'indices',[]);
tbl = nan(n,N);
for i = 1:N
    if ~isnan(st(i))
        periods(i).indices = (st(i):en(i))';
        tbl(1:numel(periods(i).indices),i) = periods(i).indices;
    end
end
tbl = tbl(any(~isnan(tbl),2),:);
if isempty(tbl), tbl = nan(1,N); end

    function pickRange(i)
        mode = questdlg(sprintf('Enter "%s" by:',names{i}),'Enter period', ...
            'Click on plot','Type values','Click on plot');
        if strcmp(mode,'Type values')
            typeRange(i); return;
        end
        while true
            title(ax,sprintf('%s: click START then END of  "%s"  (Enter = not used)',figTitle,names{i}), ...
                'Interpreter','none');
            figure(fig);
            [x,~] = ginput(2);
            if numel(x) < 2
                if strcmp(questdlg(sprintf('Mark "%s" as not used?',names{i}),'Period','Not used','Try again','Try again'),'Not used')
                    st(i) = nan; en(i) = nan; drawPatch(i); return;
                end
                continue;
            end
            x = sort(round(x));
            st(i) = max(1,x(1)); en(i) = min(n,x(2));
            if en(i) < st(i), continue; end
            drawPatch(i); return;
        end
    end

    function typeRange(i)
        if isnan(st(i)), d = {'',''}; else, d = {num2str(st(i)),num2str(en(i))}; end
        while true
            a = inputdlg({sprintf('"%s" START frame index (1-%d)  [both blank = not used]:',names{i},n), ...
                          sprintf('"%s" END frame index (1-%d):',names{i},n)}, ...
                          'Type period',1,d);
            if isempty(a)                       % cancelled: keep the current value
                return;
            end
            s1 = str2double(a{1}); e1 = str2double(a{2});
            if isempty(strtrim(a{1})) && isempty(strtrim(a{2}))
                st(i) = nan; en(i) = nan; drawPatch(i); return;
            end
            if isnan(s1) || isnan(e1)
                uiwait(msgbox('Please type numbers for both start and end.','Type period','modal')); continue;
            end
            lo = max(1,round(min(s1,e1))); hi = min(n,round(max(s1,e1)));
            st(i) = lo; en(i) = hi; drawPatch(i); return;
        end
    end

    function drawPatch(i)
        if isgraphics(hP(i)), delete(hP(i)); end
        if isgraphics(hT(i)), delete(hT(i)); end
        if isnan(st(i)), return; end
        yl = ylim(ax);
        hP(i) = patch(ax,[st(i) en(i) en(i) st(i)],[yl(1) yl(1) yl(2) yl(2)],colors(i,:), ...
            'FaceAlpha',0.15,'EdgeColor',colors(i,:),'HitTest','off');
        hT(i) = text(ax,st(i),yl(2)-(i-1)*0.04*diff(yl),[' ' names{i}],'Color',colors(i,:), ...
            'VerticalAlignment','top','FontWeight','bold','HitTest','off');
        uistack(hP(i),'bottom');
    end
end
