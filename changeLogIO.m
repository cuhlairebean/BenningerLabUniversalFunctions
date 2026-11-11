function T = changeLogIO(action, arg1, arg2)
%CHANGELOGIO  Keep a cumulative "Change log" sheet inside an analysis workbook (.xlsx).
%
%   WHY THIS EXISTS
%   Every time the analysis is run, re-run (Redo) or amended, things change: cells get removed or restored, new
%   cells are appended or all cells are redrawn, periods are re-selected, the CaMPARI threshold is redone, ...
%   This function records those changes in a sheet called 'Change log' in each output workbook that holds
%   analysis / mean values (<name>_islet#_normalized.xlsx, <name>_islet#_FreqAmpDutyCycle.xlsx,
%   Analysis_dataset_detrended.xlsx), so anyone opening the workbook can see how its numbers came to be.
%
%   THE LOG SHEET
%   One row per event, four text columns:
%       Time      when the run happened (yyyy-mm-dd HH:MM:SS)
%       File      the data set it is about (e.g. 10-29-26_5mMPC_d1_islet2)
%       Category  what kind of event (Run, Cells removed, Cells appended, Periods, Normalization, CaMPARI, ...)
%       Detail    free text (e.g. the removed cell names, the new threshold)
%   The log is CUMULATIVE: new rows are added below the rows of earlier runs; nothing is overwritten.
%
%   IMPORTANT - why there is a 'read' and a 'write'
%   Some exporters (processAndExportPeriods, freqAmpDutyCycleExport) DELETE and rewrite the whole workbook, which
%   would wipe the log sheet. So the calling script does:
%       prev = changeLogIO('read', xlsxPath);          % 1. read the old log BEFORE exporting
%       ... export the workbook (may delete it) ...
%       new  = changeLogIO('rows', 'islet label', {Category, Detail; ...});   % 2. build this run's rows
%       changeLogIO('write', xlsxPath, [prev; new]);   % 3. write old + new rows AFTER exporting
%   Only the 'Change log' sheet is written; the other sheets of the workbook are left untouched.
%
%   USAGE
%       T = changeLogIO('read',  xlsxPath)
%           Returns the existing 'Change log' as a table (0 rows if the file or sheet does not exist yet).
%       T = changeLogIO('rows',  fileLabel, C)
%           C is an N x 2 cell array {Category, Detail; ...}. Returns an N-row table stamped with the current
%           time and fileLabel, ready to be vertically concatenated with the old log.
%       changeLogIO('write', xlsxPath, T)
%           Writes table T to sheet 'Change log' of xlsxPath. Failure (e.g. the file is open in Excel) only
%           produces a warning; it never stops the analysis.
%
%   All columns are kept as MATLAB string arrays so tables read from Excel and tables built here can always be
%   concatenated.
vn = {'Time','File','Category','Detail'};
switch lower(action)
    case 'read'
        T = table(strings(0,1),strings(0,1),strings(0,1),strings(0,1),'VariableNames',vn);
        xlsxPath = arg1;
        try
            if exist(xlsxPath,'file') == 2 && any(strcmp(sheetnames(xlsxPath),'Change log'))
                R = readtable(xlsxPath,'Sheet','Change log','TextType','string');
                if all(ismember(vn,R.Properties.VariableNames))
                    R = R(:,vn);
                    for v = 1:numel(vn), R.(vn{v}) = string(R.(vn{v})); end
                    T = R;
                end
            end
        catch
            % unreadable / missing log: start a new one
        end

    case 'rows'
        fileLabel = arg1; C = arg2;
        n = size(C,1);
        stamp = string(datestr(now,'yyyy-mm-dd HH:MM:SS'));
        T = table(repmat(stamp,n,1), repmat(string(fileLabel),n,1), ...
            string(C(:,1)), string(C(:,2)), 'VariableNames',vn);

    case 'write'
        xlsxPath = arg1; T = arg2;
        try
            writetable(T,xlsxPath,'Sheet','Change log');
        catch ME
            warning('Could not write the Change log sheet of %s (%s).',xlsxPath,ME.message);
        end

    otherwise
        error('changeLogIO: unknown action ''%s'' (use ''read'', ''rows'' or ''write'').',action);
end
end
