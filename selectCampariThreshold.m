function [threshold, bg] = selectCampariThreshold(greenImg, redImg, baseName, excelPath)
% SELECTCAMPARITHRESHOLD  Pick the red/green ratio threshold separating "green" from "red"
% CaMPARI cells for ONE experiment and save it into an existing Excel file.
%
%   threshold = selectCampariThreshold(greenImg, redImg, baseName, excelPath)
%
%   greenImg, redImg : 2-D images (e.g. channels of the <name>_pre.czi)
%   baseName         : row key(s) in the Excel file: a name, or a cell array of names
%                      (the same threshold is written for each, e.g. one per islet)
%   excelPath        : existing .xlsx (e.g. CaMPARI_Color_Thresholds.xlsx)
%
% Steps (same method as campari_color_threshold_classifyer.m):
%   1) draw a BACKGROUND ROI                      (ENTER to commit)
%   2) draw GREEN cells ROIs                      (ENTER commits each, SPACE = done)
%   3) draw RED   cells ROIs                      (ENTER commits each, SPACE = done)
%   4) ratio = (red - red_bkg)/(green - green_bkg) per cell is plotted; CLICK the threshold
%   5) sheets 'Thresholds' (Filename, Threshold) and 'background_values'
%      (Filename, Green_Background, Red_Background) are updated: row for baseName
%      is replaced if it exists, otherwise appended. Other sheets are not touched.
% Window: red channel (magenta) | green channel | overlay where the ROIs are drawn (all zoom together).
% Keys: n/m = dim/brighten magenta, f/g = dim/brighten green, v/b = both; mouse wheel or =/- zoom, 0 = reset view.
threshold = NaN; bg = [NaN NaN];   % bg = [green red] background values (second output)
keys = cellstr(baseName);          % one or several Excel row keys (e.g. filename_islet1, filename_islet2)
baseName = keys{1};                % used for window titles
greenImg = double(greenImg); redImg = double(redImg);
if ~isequal(size(greenImg),size(redImg)), error('Green and red images must be the same size.'); end

fig = figure('Name',['CaMPARI threshold: ' baseName],'NumberTitle','off','Units','normalized','Position',[0.03 0.1 0.94 0.75]);
axM = axes('Parent',fig,'Position',[0.01 0.12 0.2 0.75]);    % red channel shown as magenta
axG = axes('Parent',fig,'Position',[0.22 0.12 0.2 0.75]);    % green channel
ax  = axes('Parent',fig,'Position',[0.44 0.03 0.55 0.85]);   % overlay: ROIs are drawn here
gM = 1; gG = 1;                                              % brightness (gain) of the magenta / green channels
cM = 1; cG = 1;                                              % contrast (display exponent: >1 darkens dim pixels = more contrast)
nrm = @(im) (im-min(im(:)))/max(eps,max(im(:))-min(im(:)));
rN = nrm(redImg); gN = nrm(greenImg);
dispCh = @(x,g,c) min(1,x*g).^c;
hM = imshow(cat(3,dispCh(rN,gM,cM),zeros(size(rN)),dispCh(rN,gM,cM)),'Parent',axM);
hG = imshow(cat(3,zeros(size(gN)),dispCh(gN,gG,cG),zeros(size(gN))),'Parent',axG);
hImg = imshow(cat(3,dispCh(rN,gM,cM),dispCh(gN,gG,cG),dispCh(rN,gM,cM)),'Parent',ax);
title(axM,'RED channel (magenta)','Interpreter','none'); title(axG,'GREEN channel','Interpreter','none');
linkaxes([axM axG ax]);
set(fig,'KeyPressFcn',@(~,e) keyCb(e),'WindowScrollWheelFcn',@(~,e) scrollCb(e));
    function keyCb(e)
        switch lower(e.Key)
            case 'n', gM = max(0.05,gM/1.25);          % magenta dimmer
            case 'm', gM = gM*1.25;                    % magenta brighter
            case 'f', gG = max(0.05,gG/1.25);          % green dimmer
            case 'g', gG = gG*1.25;                    % green brighter
            case 'v', gM = max(0.05,gM/1.25); gG = max(0.05,gG/1.25);   % both dimmer
            case 'b', gM = gM*1.25; gG = gG*1.25;                       % both brighter
            case 'j', cM = max(0.2,cM/1.15);           % magenta contrast down
            case 'k', cM = min(5,cM*1.15);             % magenta contrast up
            case 's', cG = max(0.2,cG/1.15);           % green contrast down
            case 'd', cG = min(5,cG*1.15);             % green contrast up
            case 'comma',  cM = max(0.2,cM/1.15); cG = max(0.2,cG/1.15);   % both contrast down
            case 'period', cM = min(5,cM*1.15);   cG = min(5,cG*1.15);     % both contrast up
            case {'equal','add'},      zoomAt(ax,[],0.8);                % zoom in
            case {'hyphen','subtract'}, zoomAt(ax,[],1.25);              % zoom out
            case '0', axis(ax,'image');                                  % reset view
            case 'r', gM = 1; gG = 1; cM = 1; cG = 1;                    % reset brightness + contrast
            otherwise, return;
        end
        refreshImgs();
    end
    function refreshImgs()
        m = dispCh(rN,gM,cM); g = dispCh(gN,gG,cG);
        set(hM,'CData',cat(3,m,zeros(size(rN)),m));
        set(hG,'CData',cat(3,zeros(size(gN)),g,zeros(size(gN))));
        set(hImg,'CData',cat(3,m,g,m));
    end
    function scrollCb(e)
        % mouse wheel zooms about the cursor in whichever panel it is over (all panels stay linked)
        for aa = [ax axG axM]
            cp = aa.CurrentPoint(1,1:2);
            xl = aa.XLim; yl = aa.YLim;
            if cp(1)>=xl(1) && cp(1)<=xl(2) && cp(2)>=yl(1) && cp(2)<=yl(2)
                zoomAt(aa,cp,1.25^sign(e.VerticalScrollCount)); return;
            end
        end
    end
    function zoomAt(aa,c,f)
        xl = aa.XLim; yl = aa.YLim;
        if isempty(c), c = [mean(xl) mean(yl)]; end
        sz = [size(redImg,2) size(redImg,1)];
        nx = c(1)+(xl-c(1))*f; ny = c(2)+(yl-c(2))*f;
        if diff(nx) >= sz(1), nx = [0.5 sz(1)+0.5]; end      % never zoom out past the full image
        if diff(ny) >= sz(2), ny = [0.5 sz(2)+0.5]; end
        aa.XLim = nx; aa.YLim = ny;
    end
keyHelp = sprintf('BRIGHTNESS: [n/m] magenta -/+ | [f/g] green -/+ | [v/b] both -/+     CONTRAST: [j/k] magenta -/+ | [s/d] green -/+ | [,/.] both -/+\nZOOM: scroll or [=]/[-] | [0] reset view | [r] reset brightness+contrast');


% 1) background
title(ax,{'BACKGROUND = AUTOFLUORESCENCE: draw an ROI on CaMPARI-NEGATIVE cells (not empty space), then press ENTER',keyHelp},'Interpreter','none');
roi = drawfreehand(ax,'Color',[1 1 0]);
waitKey(fig,{'return'});
bgMask = createMask(roi);
gBkg = mean(greenImg(bgMask)); rBkg = mean(redImg(bgMask)); bg = [gBkg rBkg];
roi.InteractionsAllowed = 'none';

% 2) + 3) cells
greenRatio = collectRatios('GREEN',[0 1 0]);
redRatio   = collectRatios('RED',[1 0 0]);
    function ratios = collectRatios(label,col)
        ratios = [];
        while isvalid(fig)
            title(ax,{sprintf('%s CELLS (#%d): draw ROI, ENTER = commit',label,numel(ratios)+1), ...
                ['SPACE = done with this group | ' keyHelp]},'Interpreter','none');
            r = drawfreehand(ax,'Color',col);
            k = waitKey(fig,{'return','space'});
            if strcmp(k,'space')
                if isvalid(r), delete(r); end
                break;
            end
            m = createMask(r);
            ratios(end+1,1) = (mean(redImg(m))-rBkg)/(mean(greenImg(m))-gBkg); %#ok<AGROW>
            r.InteractionsAllowed = 'none';
        end
    end
if isvalid(fig), close(fig); end
if isempty(greenRatio) && isempty(redRatio)
    warning('No cells were selected - nothing saved.'); return;
end

% 4) click threshold
sty = get_campari_style();
f2 = figure('Name',['Ratio threshold: ' baseName],'NumberTitle','off'); hold on;
scatter(1+(rand(size(greenRatio))-0.5)*0.25,greenRatio,60,sty.colors('green'),'filled');
scatter(2+(rand(size(redRatio))-0.5)*0.25,redRatio,60,sty.colors('magenta'),'filled');
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'Green cells','Red cells'});
ylabel('Red/Green ratio (background-subtracted)'); grid on;
fnT = sty.font_name; if ~any(strcmp(listfonts,fnT)), fnT = get(groot,'defaultAxesFontName'); end
set(gca,'FontName',fnT,'FontSize',sty.axis_size);
accepted = false;
while ~accepted
    title(gca(f2),{baseName,'CLICK on the plot to set the threshold'},'Interpreter','none');
    [~,threshold] = ginput(1);
    if isempty(threshold), threshold = NaN; break; end
    hl = yline(threshold,'b--',sprintf('Threshold = %.3f',threshold),'LineWidth',2);
    nErr = nnz(greenRatio >= threshold) + nnz(redRatio < threshold);
    nTot = numel(greenRatio)+numel(redRatio);
    a = questdlg(sprintf('Threshold = %.3f  (%d of %d cells misclassified, %.1f%%). Use it?', ...
        threshold,nErr,nTot,100*nErr/nTot),'Threshold','Yes','Re-pick','Cancel','Yes');
    if strcmp(a,'Yes'), accepted = true;
    elseif strcmp(a,'Cancel'), threshold = NaN; break;
    else, delete(hl);
    end
end
if isvalid(f2), close(f2); end
if isnan(threshold), disp('Threshold not saved.'); return; end

% 5) save into the existing Excel file
for q = 1:numel(keys)
    upsertRow(excelPath,'Thresholds',{'Filename','Threshold'},keys{q},threshold);
    upsertRow(excelPath,'background_values',{'Filename','Green_Background','Red_Background'}, ...
        keys{q},[gBkg rBkg]);
    fprintf('Threshold %.3f for %s saved to %s\n',threshold,keys{q},excelPath);
end
end

function k = waitKey(fig,keys)
% wait until one of the keys is pressed in fig (a/s are handled by the figure callback)
k = '';
old = fig.KeyPressFcn;
state.k = '';
set(fig,'KeyPressFcn',@(s,e) cb(s,e));
    function cb(s,e)
        if ~isempty(old), old(s,e); end
        if ismember(lower(e.Key),keys), state.k = lower(e.Key); uiresume(s); end
    end
while isvalid(fig) && isempty(state.k)
    uiwait(fig);
end
k = state.k;
if isvalid(fig), set(fig,'KeyPressFcn',old); end
end

function upsertRow(excelPath,sheet,varNames,baseName,values)
% replace the row whose Filename == baseName, or append it, then rewrite the sheet
T = table();
if exist(excelPath,'file')
    try
        T = readtable(excelPath,'Sheet',sheet);
        T.Properties.VariableNames = varNames(1:width(T));
    catch
        T = table();
    end
end
newRow = [{baseName}, num2cell(values)];
if isempty(T) || height(T) == 0
    T = cell2table(newRow,'VariableNames',varNames);
else
    names = cellstr(string(T.(varNames{1})));
    i = find(strcmp(names,baseName),1);
    if isempty(i), i = height(T)+1; end
    T(i,:) = cell2table(newRow,'VariableNames',varNames);
end
writetable(T,excelPath,'Sheet',sheet);
end
