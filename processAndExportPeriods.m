function written = processAndExportPeriods(CellTC, periods, timeVec, outPrefix, procChoices, cellLabels)
%PROCESSANDEXPORTPERIODS  Detrend / baseline-normalize cell timecourses and export to Excel.
%
%   written = processAndExportPeriods(CellTC, periods, timeVec, outPrefix, procChoices, cellLabels)
%
%   CellTC       time x cells matrix (one islet)
%   periods      struct array with .name .start .end  (1-based frame indices, as saved by
%                selectTimeRanges / *_PeriodIndices.mat). 'Baseline' is needed for normalization;
%                'Low glucose','High glucose','Drug' only LABEL frames and give per-period means.
%   timeVec      time stamp of every frame (seconds) - frame index is used if empty
%   outPrefix    e.g. fullfile(output_dir,[filename '_islet2'])  ->  <outPrefix>_<processing>.xlsx
%   procChoices  cellstr subset of {'detrended','normalized','normalized_inverted',
%                'normalized_minus_mean'}   (default: all)
%   cellLabels   cellstr, one name per column of CellTC (default cell1..cellN)
%
%   The WHOLE matrix is processed (every frame, including frames outside any selected period):
%     detrended              one linear trendline per cell over the whole recording removed
%     normalized             (x - mean_baseline) / mean_baseline   (ORIGINAL data)
%     normalized_inverted    -1 * normalized
%     normalized_minus_mean  normalized - mean(normalized over the whole recording)
%   mean_baseline = mean of each cell over the 'Baseline' period you selected.
%
%   Excel sheets:  'data'  index | time | condition | one column per cell, all frames. The
%                          'condition' column only labels which selected period(s) a frame is in
%                          (blank = in none).
%                  'Means' mean of each cell over the whole recording and in each selected period
%                  'Conditions' start/end frames of the selected periods (incl. baseline)

exportNames = {'Low glucose','High glucose','Drug'};
allProcs = {'detrended','normalized','normalized_inverted','normalized_minus_mean'};
if nargin < 5 || isempty(procChoices), procChoices = allProcs; end
[nT, nC] = size(CellTC);
if nargin < 6 || isempty(cellLabels)
    cellLabels = arrayfun(@(j) sprintf('cell%d',j),1:nC,'UniformOutput',false);
end
if nargin < 3 || isempty(timeVec), timeVec = (1:nT)'; end
timeVec = timeVec(:);
if numel(timeVec) < nT, timeVec(end+1:nT,1) = NaN; end
timeVec = timeVec(1:nT);
written = {};

% ---- resolve ranges (clip / skip anything out of range) ----
exportR = struct('name',{},'start',{},'end',{},'note',{});
for i = 1:numel(exportNames)
    r = getRange(periods, exportNames{i}, nT);
    if ~isempty(r), exportR(end+1) = r; end %#ok<AGROW>
end
baseR = getRange(periods,'Baseline',nT);
if isempty(exportR)
    warning('processAndExportPeriods:noPeriods','No Low/High glucose or Drug period selected - the whole matrix is still processed.');
end

% ---- condition label per frame ----
condCol = repmat({''},nT,1);
for i = 1:numel(exportR)
    for f = exportR(i).start:exportR(i).end
        if isempty(condCol{f}), condCol{f} = exportR(i).name;
        else, condCol{f} = [condCol{f} ';' exportR(i).name]; end
    end
end

% ---- baseline mean per cell ----
baseMean = nan(1,nC);
if ~isempty(baseR)
    baseMean = mean(CellTC(baseR.start:baseR.end,:),1,'omitnan');
end

% ---- Conditions sheet (periods are labels only; processing uses the whole matrix) ----
rows = [exportR, baseR];
condNames = {exportR.name}';
if ~isempty(baseR), condNames{end+1,1} = 'BASELINE (normalization only)'; end
if isempty(rows)
    condTbl = table('Size',[0 7],'VariableTypes',{'cell','double','double','double','double','double','cell'}, ...
        'VariableNames',{'Period','StartIndex','EndIndex','NFrames','StartTime','EndTime','Note'});
else
    sIdx = [rows.start]'; eIdx = [rows.end]';
    condTbl = table(condNames, sIdx, eIdx, eIdx-sIdx+1, timeVec(sIdx), timeVec(eIdx), {rows.note}', ...
        'VariableNames',{'Period','StartIndex','EndIndex','NFrames','StartTime','EndTime','Note'});
end

tAll = timeVec;
if all(isnan(tAll)), tAll = (1:nT)'; end

for p = 1:numel(procChoices)
    proc = procChoices{p};
    if ~strcmp(proc,'detrended') && all(isnan(baseMean))
        warning('processAndExportPeriods:noBaseline', ...
            '%s: no usable Baseline period - skipping %s.',outPrefix,proc);
        continue
    end
    % ---- whole matrix, every frame (including frames outside any selected period) ----
    switch proc
        case 'detrended'
            M = nan(nT,nC);
            for c = 1:nC
                y = CellTC(:,c); ok = ~isnan(y) & ~isnan(tAll);
                if nnz(ok) >= 2
                    co = polyfit(tAll(ok),y(ok),1);      % one linear trend per cell over the whole recording
                    M(:,c) = y - polyval(co,tAll);
                else
                    M(:,c) = y;
                end
            end
        otherwise
            M = (CellTC - baseMean) ./ baseMean;         % ORIGINAL data / baseline mean
            M(:, baseMean == 0 | isnan(baseMean)) = NaN;
            if strcmp(proc,'normalized_inverted')
                M = -1 * M;
            elseif strcmp(proc,'normalized_minus_mean')
                M = M - mean(M,1,'omitnan');             % mean of the whole recording
            end
    end

    % ---- Means: whole recording first, then each selected period ----
    meanNames = [{'Whole recording'}, {exportR.name}];
    meansM = nan(numel(meanNames),nC);
    meansM(1,:) = mean(M,1,'omitnan');
    for i = 1:numel(exportR)
        meansM(i+1,:) = mean(M(exportR(i).start:exportR(i).end,:),1,'omitnan');
    end

    dataTbl = table((1:nT)',timeVec,condCol,'VariableNames',{'index','time','condition'});
    dataTbl = [dataTbl, array2table(M,'VariableNames',cellLabels)];
    meansTbl = [table(meanNames(:),'VariableNames',{'Period'}), ...
                array2table(meansM,'VariableNames',cellLabels), ...
                table(mean(meansM,2,'omitnan'),'VariableNames',{'mean_across_cells'})];

    xlsxPath = [outPrefix '_' proc '.xlsx'];
    if exist(xlsxPath,'file'), delete(xlsxPath); end
    writetable(dataTbl,  xlsxPath,'Sheet','data');
    writetable(meansTbl, xlsxPath,'Sheet','Means');
    writetable(condTbl,  xlsxPath,'Sheet','Conditions');
    written{end+1} = xlsxPath; %#ok<AGROW>
    fprintf('  saved %s\n',xlsxPath);
end

end

function r = getRange(periods, name, nT)
%GETRANGE  Valid [start end] for a named period, or [] if missing/unusable.
r = [];
j = find(strcmpi({periods.name},name),1);
if isempty(j), return; end
s = periods(j).start; e = periods(j).end;
if isempty(s) || isempty(e) || isnan(s) || isnan(e), return; end
note = '';
if s > e, warning('%s: start %d > end %d - skipped.',name,s,e); return; end
if s < 1, note = [note sprintf('start %d clipped to 1; ',s)]; s = 1; end
if s > nT, warning('%s: start %d is beyond the %d frames - skipped.',name,s,nT); return; end
if e > nT, note = [note sprintf('end %d clipped to %d; ',e,nT)]; e = nT; end
r = struct('name',name,'start',s,'end',e,'note',strtrim(note));
end
