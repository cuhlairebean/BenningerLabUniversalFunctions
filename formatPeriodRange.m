function [idxStr, timeStr, startIdx, startTime] = formatPeriodRange(periods, name, timeVec, nFrames)
% FORMATPERIODRANGE  Text for one selected period, as used in the analysis dataset sheet.
%   [idxStr, timeStr, startIdx, startTime] = formatPeriodRange(periods, name, timeVec, nFrames)
%   periods  struct array (.name .start .end) from selectTimeRanges
%   timeVec  time stamp (s) of every frame (frame index is used if empty)
%   nFrames  number of frames; an end index >= nFrames is written as 'end'
%   idxStr   e.g. '750-1250' or '2200-end'      (blank if the period is not used)
%   timeStr  e.g. '761.31-1268.17'  (seconds, 2 decimals)
%   startIdx / startTime : first index and its time (NaN if not used)
idxStr = ''; timeStr = ''; startIdx = NaN; startTime = NaN;
j = find(strcmpi({periods.name},name),1);
if isempty(j), return; end
s = periods(j).start; e = periods(j).end;
if isempty(s) || isempty(e) || isnan(s) || isnan(e), return; end
timeVec = timeVec(:);
if isempty(timeVec), timeVec = (1:nFrames)'; end
tAt = @(i) timeVec(min(max(i,1),numel(timeVec)));
startIdx = s; startTime = tAt(s);
if e >= nFrames, endTxt = 'end'; e = nFrames; else, endTxt = sprintf('%d',e); end
idxStr  = sprintf('%d-%s',s,endTxt);
timeStr = sprintf('%.2f-%.2f',tAt(s),tAt(e));
end
