function initial = loadPreviousPeriods(file, periodNames)
%LOADPREVIOUSPERIODS  Read earlier period selections (start/end frame per period).
%   initial = loadPreviousPeriods(file, periodNames)
%   file may be  *_PeriodIndices.mat  (IsletPeriods struct)  or  *_PeriodRanges.csv.
%   Returns an Nx2 matrix [start end] in the order of periodNames (NaN = not used).
N = numel(periodNames);
initial = nan(N,2);
[~,~,ext] = fileparts(file);
if strcmpi(ext,'.mat')
    S = load(file,'IsletPeriods');
    names = {S.IsletPeriods.name};
    st = [S.IsletPeriods.start]; en = [S.IsletPeriods.end];
else
    T = readtable(file);
    names = cellstr(string(T{:,1}));
    st = T{:,2}'; en = T{:,3}';
end
for i = 1:N
    j = find(strcmpi(strtrim(names),periodNames{i}),1);
    if ~isempty(j), initial(i,:) = [st(j) en(j)]; end
end
end
