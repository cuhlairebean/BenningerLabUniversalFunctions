function threshold = selectCampariThreshold(greenImg, redImg, baseName, excelPath)
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
% Image keys: a = brighten, s = dim.
threshold = NaN;
keys = cellstr(baseName);          % one or several Excel row keys (e.g. filename_islet1, filename_islet2)
baseName = keys{1};                % used for window titles
greenImg = double(greenImg); redImg = double(redImg);
if ~isequal(size(greenImg),size(redImg)), error('Green and red images must be the same size.'); end

fig = figure('Name',['CaMPARI threshold: ' baseName],'NumberTitle','off');
ax = axes('Parent',fig);
gain = 1;
rgbBase = cat(3, redImg/max(redImg(:)), greenImg/max(greenImg(:)), zeros(size(redImg)));
hImg = imshow(min(1,rgbBase*gain),'Parent',ax);
set(fig,'KeyPressFcn',@(~,e) keyCb(e));
    function keyCb(e)
        switch lower(e.Key)
            case 'a', gain = gain*1.25; set(hImg,'CData',min(1,rgbBase*gain));
            case 's', gain = max(0.05,gain/1.25); set(hImg,'CData',min(1,rgbBase*gain));
        end
    end

% 1) background
title(ax,{'BACKGROUND: draw an ROI, then press ENTER','[a] brighten | [s] dim'},'Interpreter','none');
roi = drawfreehand(ax,'Color',[1 1 0]);
waitKey(fig,{'return'});
bgMask = createMask(roi);
gBkg = mean(greenImg(bgMask)); rBkg = mean(redImg(bgMask));
roi.InteractionsAllowed = 'none';

% 2) + 3) cells
greenRatio = collectRatios('GREEN',[0 1 0]);
redRatio   = collectRatios('RED',[1 0 0]);
    function ratios = collectRatios(label,col)
        ratios = [];
        while isvalid(fig)
            title(ax,{sprintf('%s CELLS (#%d): draw ROI, ENTER = commit',label,numel(ratios)+1), ...
                'SPACE = done with this group | [a] brighten | [s] dim'},'Interpreter','none');
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
f2 = figure('Name',['Ratio threshold: ' baseName],'NumberTitle','off'); hold on;
scatter(1+(rand(size(greenRatio))-0.5)*0.25,greenRatio,60,[0 0.8 0],'filled');
scatter(2+(rand(size(redRatio))-0.5)*0.25,redRatio,60,[0.9 0 0],'filled');
xlim([0.5 2.5]); xticks([1 2]); xticklabels({'Green cells','Red cells'});
ylabel('Red/Green ratio (background-subtracted)'); grid on;
accepted = false;
while ~accepted
    title(f2,{baseName,'CLICK on the plot to set the threshold'},'Interpreter','none');
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
