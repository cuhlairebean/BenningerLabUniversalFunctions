function DataOut = freqAmpDutyCycleExport(NormData, periods, timeVec, cellLabels, xlsxPath)
%FREQAMPDUTYCYCLEEXPORT  Frequency, amplitude and duty cycle per cell for 4 phases -> Excel.
%
%   DataOut = freqAmpDutyCycleExport(NormData, periods, timeVec, cellLabels, xlsxPath)
%
%   NormData    time x cells normalized matrix (e.g. mats.normalized from processAndExportPeriods)
%   periods     struct array (.name .start .end) as saved by selectTimeRanges
%   timeVec     time stamp (s) of every frame; Fs = 1/median(diff(timeVec)). Empty -> Fs = 1
%   cellLabels  cellstr, one per column
%   xlsxPath    output .xlsx
%
%   Phases (period names) -> column headers:
%       'Low glucose' -> '2mM', 'First phase', 'Second phase', 'Drug'
%   Workbook layout follows 'Example Output Excel .xlsx':
%       sheets Frequency | Amplitude | Duty Cycle (rows = cells, columns = phases) and Means.
%   Math (smoothing, FFT frequency/amplitude, duty cycle) is the lab's DUTY CYCLE / FREQUENCY /
%   AMPLITUDE script.

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

DataOut.dutycycle = nan(nC,nPh);
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
        normData = normalize(detrend(y));
        DataOut.dutycycle(c,p) = sum(normData >= thresh)/length(normData);

        sg = detrend(y);
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

% ---- Excel (layout like the example workbook) ----
if exist(xlsxPath,'file'), delete(xlsxPath); end
sheets = {'Frequency',DataOut.freqvec; 'Amplitude',DataOut.amp; 'Duty Cycle',DataOut.dutycycle};
for k = 1:size(sheets,1)
    writecell([headers; num2cell(sheets{k,2})],xlsxPath,'Sheet',sheets{k,1});
    % cell names in the column after the data
    writecell([{'Cell'}; cellLabels(:)],xlsxPath,'Sheet',sheets{k,1},'Range',sprintf('%s1',colLetter(nPh+2)));
end
meansCell = [{''},headers; ...
    {'Frequency'},  num2cell(mean(DataOut.freqvec,1,'omitnan')); ...
    {'Amplitude'},  num2cell(mean(DataOut.amp,1,'omitnan')); ...
    {'Duty Cycle'}, num2cell(mean(DataOut.dutycycle,1,'omitnan'))];
writecell(meansCell,xlsxPath,'Sheet','Means');
fprintf('  saved %s\n',xlsxPath);
end

function s = colLetter(n)
s = '';
while n > 0
    r = mod(n-1,26);
    s = [char(65+r) s]; %#ok<AGROW>
    n = floor((n-1)/26);
end
end
