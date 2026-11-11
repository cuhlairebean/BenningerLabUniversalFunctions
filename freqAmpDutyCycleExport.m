function DataOut = freqAmpDutyCycleExport(NormData, periods, timeVec, cellLabels, xlsxPath, ratio, cls)
%FREQAMPDUTYCYCLEEXPORT  Frequency, amplitude and duty cycle per cell for 4 phases -> Excel.
%
%   DataOut = freqAmpDutyCycleExport(NormData, periods, timeVec, cellLabels, xlsxPath, ratio, cls)
%
%   NormData    time x cells normalized matrix (e.g. mats.normalized from processAndExportPeriods)
%   periods     struct array (.name .start .end) as saved by selectTimeRanges
%   timeVec     time stamp (s) of every frame; Fs = 1/median(diff(timeVec)). Empty -> Fs = 1
%   cellLabels  cellstr, one per column
%   xlsxPath    output .xlsx
%   ratio, cls  OPTIONAL (CaMPARI): red:green ratio per cell (1 x cells) and its class cellstr ('red'/'green'/
%               'unclassified'). They are written as two extra columns on every sheet, and the Means sheet gets
%               extra rows with the mean of the green cells and of the red cells for every measure.
%
%   Phases (period names) -> column headers:
%       'Low glucose' -> '2mM', 'First phase', 'Second phase', 'Drug'
%   Workbook layout follows 'Example Output Excel .xlsx':
%       sheets Frequency | Amplitude | Duty Cycle (rows = cells, columns = phases) and Means.
%   Math (smoothing, FFT frequency/amplitude, duty cycle) is the lab's DUTY CYCLE / FREQUENCY /
%   AMPLITUDE script.

if nargin < 6, ratio = []; end
if nargin < 7, cls = {}; end
periodNames = {'Low glucose','First phase','Second phase','Drug'};
headers     = {'2mM','First Phase','Second Phase','Drug'};
nPh = numel(periodNames);
nC  = size(NormData,2);
nT  = size(NormData,1);

timeVec = timeVec(:);
if numel(timeVec) >= 2 && ~all(isnan(timeVec))
    Fs = 1/median(diff(timeVec),'omitnan');
else
    Fs = 1;
    warning('freqAmpDutyCycleExport:noTime','No time stamps - Fs = 1 (frequency is in cycles/frame).');
end
fprintf('  Frequency analysis: Fs = %.4f Hz\n',Fs);

neighbors  = 3;                 % moving-average smoothing
windowSize = 2*neighbors + 1;
thresh     = 0.3;               % duty cycle "on" threshold (z-score)

% DUTY CYCLE - same two methods as the Python pipeline (NormalizeandVisualize):
%   M1 = fraction of frames at/above the (detrended) mean              -> used for activity at HIGH glucose / the other phases
%   M2 = fraction of frames at/above mean + 2 standard deviations      -> used for BASELINE (2mM), only significant events count
% Both are computed on the smoothed (7-frame moving average) + linearly detrended trace of each phase.
% The 'Duty Cycle' sheet / DataOut.dutycycle uses M2 for the 2mM (baseline) column and M1 for First Phase, Second Phase
% and Drug. The full M1 and M2 results for every phase are also written ('... (all phases)' sheets).
DataOut.dutycycleM1 = nan(nC,nPh);
DataOut.dutycycleM2 = nan(nC,nPh);
DataOut.freqvec   = nan(nC,nPh);
DataOut.amp       = nan(nC,nPh);

for p = 1:nPh
    j = find(strcmpi({periods.name},periodNames{p}),1);
    if isempty(j) || isempty(periods(j).start) || isempty(periods(j).end) ...
            || isnan(periods(j).start) || isnan(periods(j).end)
        warning('Period "%s" not selected - left blank.',periodNames{p});
        continue
    end
    s = max(1,periods(j).start); e = min(nT,periods(j).end);
    if e - s < 3, warning('Period "%s" too short - left blank.',periodNames{p}); continue; end
    X = movmean(NormData(s:e,:),windowSize,1);
    for c = 1:nC
        y = X(:,c);
        if any(isnan(y)), continue; end
        sg = detrend(y);
        DataOut.dutycycleM1(c,p) = dutyCycleM1(sg);
        DataOut.dutycycleM2(c,p) = dutyCycleM2(sg);
        L = length(sg);
        fts = fft(sg - mean(sg));
        nHalf = floor(L/2) + 1;
        f = Fs*(0:nHalf-1)/L;
        P1 = abs(fts/L);
        P1 = P1(1:nHalf);
        P1(2:end-1) = 2*P1(2:end-1);
        [~,idx] = max(P1);
        DataOut.freqvec(c,p) = f(idx);
        DataOut.amp(c,p)     = abs(fts(idx));
    end
end
% combined duty cycle: M2 for the baseline (2mM = phase 1), M1 for all other phases
DataOut.dutycycle = DataOut.dutycycleM1;
DataOut.dutycycle(:,1) = DataOut.dutycycleM2(:,1);

% ---- Excel (layout like the example workbook) ----
% CaMPARI (optional ratio / cls): every sheet also gets the cell's red:green ratio and its class (red / green /
% unclassified) in the two columns right after the 'Cell' column, and the Means sheet gets the mean of the
% GREEN cells and of the RED cells next to the mean of ALL cells.
hasCls = ~isempty(ratio) && ~isempty(cls) && numel(ratio) == nC && numel(cls) == nC;
if exist(xlsxPath,'file'), delete(xlsxPath); end
sheets = {'Frequency',DataOut.freqvec; 'Amplitude',DataOut.amp; 'Duty Cycle',DataOut.dutycycle; ...
    'Duty Cycle M1 (all phases)',DataOut.dutycycleM1; 'Duty Cycle M2 (all phases)',DataOut.dutycycleM2};
for k = 1:size(sheets,1)
    writecell([headers; num2cell(sheets{k,2})],xlsxPath,'Sheet',sheets{k,1});
    % cell names in the column after the data
    writecell([{'Cell'}; cellLabels(:)],xlsxPath,'Sheet',sheets{k,1},'Range',sprintf('%s1',colLetter(nPh+2)));
    if hasCls
        writecell([{'RedGreenRatio','Class'}; num2cell(ratio(:)), cls(:)],xlsxPath,'Sheet',sheets{k,1}, ...
            'Range',sprintf('%s1',colLetter(nPh+3)));
    end
end
metrics = sheets;
meansCell = [{''},headers];
for m = 1:size(metrics,1)
    meansCell = [meansCell; {metrics{m,1}}, num2cell(mean(metrics{m,2},1,'omitnan'))]; %#ok<AGROW>
    if hasCls
        isG = strcmp(cls(:),'green'); isR = strcmp(cls(:),'red');
        meansCell = [meansCell; ...
            {[metrics{m,1} ' (green cells, n=' num2str(nnz(isG)) ')']}, num2cell(mean(metrics{m,2}(isG,:),1,'omitnan')); ...
            {[metrics{m,1} ' (red cells, n=' num2str(nnz(isR)) ')']},   num2cell(mean(metrics{m,2}(isR,:),1,'omitnan'))]; %#ok<AGROW>
    end
end
writecell(meansCell,xlsxPath,'Sheet','Means');
fprintf('  saved %s\n',xlsxPath);
end

function d = dutyCycleM1(series)
% Method 1: % of time the signal is above its (detrended) mean. Python: series/np.std(series); on = mean(norm).
% (dividing by the SD does not change the result, it is kept to mirror the Python code)
series = series(:);
if numel(series) < 2 || std(series,1) == 0, d = NaN; return; end
nrm = series / std(series,1);
d = sum(nrm >= mean(nrm)) / numel(series);
end

function d = dutyCycleM2(series)
% Method 2: % of time the signal is above mean + 2 standard deviations (baseline: only significant events).
% std(...,1) = population SD like numpy's np.std.
series = series(:);
if numel(series) < 2 || std(series,1) == 0, d = NaN; return; end
nrm = series / std(series,1);
d = sum(nrm >= mean(nrm) + 2*std(nrm,1)) / numel(series);
end

function s = colLetter(n)
s = '';
while n > 0
    r = mod(n-1,26);
    s = [char(65+r) s]; %#ok<AGROW>
    n = floor((n-1)/26);
end
end
