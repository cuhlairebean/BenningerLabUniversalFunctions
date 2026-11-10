function [written, mats] = processAndExportPeriods(CellTC, periods, timeVec, outPrefix, cellLabels, multiplyNeg1, writeFiles)
%PROCESSANDEXPORTPERIODS  Per-cell detrend, baseline-normalize the original data, optional x -1,
%   and export to <outPrefix>_normalized.xlsx.
%
%   [written, mats] = processAndExportPeriods(CellTC, periods, timeVec, outPrefix, cellLabels, multiplyNeg1, writeFiles)
%
%   CellTC       time x cells matrix (one islet), RAW intensities (NOT already multiplied by -1)
%   periods      struct array with .name .start .end  (1-based frame indices, as saved by
%                selectTimeRanges / *_PeriodIndices.mat). 'Baseline' is needed for normalization;
%                'Low glucose','High glucose','Drug' only LABEL frames and give per-period means.
%   timeVec      time stamp of every frame (seconds) - frame index is used if empty
%   outPrefix    e.g. fullfile(output_dir,[filename '_islet2'])  ->  <outPrefix>_normalized.xlsx
%   cellLabels   cellstr, one name per column of CellTC (default cell1..cellN)
%   multiplyNeg1 (default false) multiply the normalized data by -1 (done LAST)
%   writeFiles   (default true) false = only compute: no Excel file is written
%
%   Order of operations, on the WHOLE matrix (every frame, also outside any selected period):
%     1) detrend    PER CELL: each cell's own linear regression vs time gives
%                   trendline = slope*time + intercept, subtracted from that cell  -> .detrended
%     2) normalize  normalized = (D - b) / F0, per cell, where
%                     D  = the detrended trace,
%                     b  = mean of D over the user-selected 'Baseline' period,
%                     F0 = mean of the cell's RAW trace over the same 'Baseline' period (its own brightness)
%     3) x -1       only if multiplyNeg1, applied LAST to the normalized data         -> .normalized
%   (No condition-by-condition segmentation as in the old Python code.)
%
%   WHY THESE STEPS
%     Detrend:   every cell holds a different amount of indicator and sits at a slightly different depth, so
%                each one bleaches/drifts at its own rate. Subtracting each cell's OWN straight-line trend
%                removes that slow drift so that only the calcium activity remains.
%     (D - b):   shifts the trace so the baseline sits at exactly 0; changes are then "above / below baseline".
%     / F0:      converts the change into a FRACTION of the cell's own resting brightness (dF/F), so bright and
%                dim cells, or cells at different depths, can be compared and averaged fairly. F0 is taken from
%                the RAW trace on purpose: after detrending the cell's brightness has been subtracted away,
%                so the detrended baseline mean is only the baseline's offset from the fitted line (it can be
%                ~0 or negative for a dim cell that does not bleach, which would amplify or flip that cell).
%                Set NORM_DENOM = 'detrended' below to restore that older behaviour.
%     x -1:      optional, for signals that DECREASE when calcium rises (e.g. Fura-type indicators) so that
%                calcium rises are plotted upward. Done last so it only flips the sign, never the scaling.
%
%   mats         struct: .detrended, .trendline, .normalized (time x cells)
%   Excel sheets:  'data'      index | time | condition | one column per cell (normalized, final)
%                  'detrended' index | time | one column per cell (step 1)
%                  'Means'     mean of each cell over the whole recording and in each selected period
%                  'Conditions' start/end frames of the selected periods (incl. baseline)

exportNames = {'Low glucose','High glucose','Drug'};
[nT, nC] = size(CellTC);
if nargin < 5 || isempty(cellLabels)
    cellLabels = arrayfun(@(j) sprintf('cell%d',j),1:nC,'UniformOutput',false);
end
if nargin < 6 || isempty(multiplyNeg1), multiplyNeg1 = false; end
if nargin < 7 || isempty(writeFiles), writeFiles = true; end
if nargin < 3 || isempty(timeVec), timeVec = (1:nT)'; end
timeVec = timeVec(:);
if numel(timeVec) < nT, timeVec(end+1:nT,1) = NaN; end
timeVec = timeVec(1:nT);
written = {};
mats = struct();

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

% ---- 1) DETREND, PER CELL (same method as the Python pipeline) ----
% Each cell contains a unique concentration of the fluorescent indicator and is positioned at a slightly
% different depth, leading to varying levels of excitation light exposure and unique bleaching kinetics.
% Detrending on a per-cell basis captures these distinct baseline decays.
%   slope/intercept = linear regression of the cell's own trace vs time   (linregress)
%   trendline       = slope * time + intercept                            (y = m x + b)
%   detrended       = original array - trendline   (drift_corrected_array in Python)
D = nan(nT,nC);
trendline = nan(nT,nC);
for c = 1:nC
    y = CellTC(:,c);
    ok = ~isnan(y) & ~isnan(tAll);
    if nnz(ok) >= 2
        co = polyfit(tAll(ok),y(ok),1);              % co(1) = slope, co(2) = intercept
        tl = co(1)*tAll + co(2);
        trendline(:,c) = tl;
        D(:,c) = y - tl;
    else
        warning('processAndExportPeriods:noTrend','%s %s: not enough valid points for a trendline - left undetrended.',outPrefix,cellLabels{c});
        trendline(:,c) = 0;
        D(:,c) = y;
    end
end
mats.detrended = D;
mats.trendline = trendline;

% ---- 2) NORMALIZE the DETRENDED data to the mean of its (detrended) baseline period ----
% normalized = (D - baseline_mean) / baseline_mean, D = detrended data, baseline_mean = mean of D over the
% user-selected 'Baseline' period. (The Python condition-by-condition segmentation is intentionally not
% used: it was a limitation of the old code; the baseline is the period selected in MATLAB.)
if isempty(baseR)
    warning('processAndExportPeriods:noBaseline', ...
        '%s: no usable Baseline period - cannot normalize.',outPrefix);
    return
end
baseMean = mean(D(baseR.start:baseR.end,:),1,'omitnan');   % detrended baseline offset (numerator reference)
% ---- DENOMINATOR: each cell is normalized to ITS OWN baseline brightness F0 ----
% F0 = mean of that cell's raw trace over the selected Baseline period (always > 0 for fluorescence).
% (Dividing by the baseline mean of the DETRENDED trace, as before, is NOT a brightness: detrending removes the
%  cell's level, so that number is just the baseline's offset from the fitted line - tiny/negative for a dim
%  cell that does not bleach, which amplified or flipped it.)
% Set NORM_DENOM = 'detrended' to get the old behaviour back.
NORM_DENOM = 'rawBaseline';
rawBase = mean(CellTC(baseR.start:baseR.end,:),1,'omitnan');
switch NORM_DENOM
    case 'rawBaseline', denom = rawBase;
    otherwise,          denom = baseMean;
end
% Denominator guard against 0 / tiny values: 0.001 is added (in the direction of its sign).
DENOM_EPS = 0.001;
tiny = abs(denom) < DENOM_EPS;                        % includes exactly 0
sgn = ones(size(denom)); sgn(denom < 0) = -1;
denom(tiny) = denom(tiny) + sgn(tiny)*DENOM_EPS;
if writeFiles
    fprintf('\n[%s] normalization denominators (%s):\n',outPrefix,NORM_DENOM);
    fprintf('   %-10s  %14s  %16s\n','cell','F0 (raw base)','detrended base');
    for cc = 1:nC
        fprintf('   %-10s  %14.4g  %16.4g\n',cellLabels{cc},rawBase(cc),baseMean(cc));
    end
    fprintf('\n');
end
if any(tiny)
    warning('processAndExportPeriods:denomCorrected','%s: 0.001 added to the baseline denominator of %d cell(s): %s', ...
        outPrefix,nnz(tiny),strjoin(cellLabels(tiny),', '));
end
M = (D - baseMean) ./ denom;                          % numerator: detrended data minus its own baseline mean (=0 at baseline)
M(:, isnan(baseMean) | isnan(denom)) = NaN;

% ---- 3) multiply by -1 LAST, if selected ----
if multiplyNeg1, M = -1 * M; end
mats.normalized = M;
if ~writeFiles, return; end

% ---- Means: whole recording first, then each selected period ----
meanNames = [{'Whole recording'}, {exportR.name}];
meansM = nan(numel(meanNames),nC);
meansM(1,:) = mean(M,1,'omitnan');
for i = 1:numel(exportR)
    meansM(i+1,:) = mean(M(exportR(i).start:exportR(i).end,:),1,'omitnan');
end

dataTbl = table((1:nT)',timeVec,condCol,'VariableNames',{'index','time','condition'});
dataTbl = [dataTbl, array2table(M,'VariableNames',cellLabels)];
detTbl  = [table((1:nT)',timeVec,'VariableNames',{'index','time'}), array2table(D,'VariableNames',cellLabels)];
trendTbl = [table((1:nT)',timeVec,'VariableNames',{'index','time'}), array2table(trendline,'VariableNames',cellLabels)];
meansTbl = [table(meanNames(:),'VariableNames',{'Period'}), ...
            array2table(meansM,'VariableNames',cellLabels), ...
            table(mean(meansM,2,'omitnan'),'VariableNames',{'mean_across_cells'})];

xlsxPath = [outPrefix '_normalized.xlsx'];
if exist(xlsxPath,'file'), delete(xlsxPath); end
writetable(dataTbl,  xlsxPath,'Sheet','data');
writetable(detTbl,   xlsxPath,'Sheet','detrended');
writetable(trendTbl, xlsxPath,'Sheet','trendline');
writetable(meansTbl, xlsxPath,'Sheet','Means');
writetable(condTbl,  xlsxPath,'Sheet','Conditions');
written{end+1} = xlsxPath;
fprintf('  saved %s\n',xlsxPath);
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
