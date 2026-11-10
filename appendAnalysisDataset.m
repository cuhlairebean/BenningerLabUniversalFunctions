function appendAnalysisDataset(xlsxPath, rowKey, summary, dataTbl)
% APPENDANALYSISDATASET  Add/replace one row in sheet 1 of the analysis workbook and write the
% matrices of that row's key to its own sheet.
%
%   appendAnalysisDataset(xlsxPath, rowKey, summary, dataTbl)
%
%   xlsxPath : workbook (created if it does not exist, e.g. Analysis_dataset_detrended.xlsx)
%   rowKey   : text in the 'Filename' column AND the sheet name (e.g. '10-3-25-2mMF4_islet1')
%   summary  : struct whose field names are sheet-1 column headers (spaces written as '_' in the
%              field name are matched to either 'Glucose Conditions' or 'Glucose_Conditions').
%              Fields not matching a header are appended as new columns on the right.
%   dataTbl  : table written to the sheet named rowKey (index, time, detrended + normalized cells)
%
% Sheet 1 is read and written as raw cells, so columns/values you already have are kept.
% A row with the same Filename is replaced, otherwise a new row is appended.

defaultHeader = {'Filename','Islet_1','Islet_2','Islet_3','Islet_4','Islet_5','Islet_6','Islet_7','Islet_8', ...
    'BL_idx','HG_idx','FP_idx','SP_idx','DRG_idx','DRG_ID','Glucose Conditions','Sex','Age', ...
    'BL','HG','FP','SP','DRG'};

% ---------- sheet 1 ----------
if exist(xlsxPath,'file')
    sheets = sheetnames(xlsxPath);
    C = readcell(xlsxPath,'Sheet',sheets(1));
    if isempty(C), C = defaultHeader; end
else
    sheets = strings(0,1);
    C = defaultHeader;
end
C = cellfun(@clean,C,'UniformOutput',false);       % missing -> []
header = cellfun(@(x) char(string(x)),C(1,:),'UniformOutput',false);
norm = @(x) lower(regexprep(x,'[\s_]+',''));

fn = fieldnames(summary);
keyCol = find(strcmp(norm(header),'filename'),1);
if isempty(keyCol), header = [{'Filename'} header]; C = [cell(size(C,1),1) C]; C(1,:) = header; keyCol = 1; end
colOf = zeros(1,numel(fn));
for i = 1:numel(fn)
    j = find(strcmp(norm(header),norm(fn{i})),1);
    if isempty(j)                                    % new column on the right
        header{end+1} = strrep(fn{i},'_',' '); %#ok<AGROW>
        if startsWith(fn{i},'Islet') || contains(fn{i},'_idx') || strcmp(fn{i},'DRG_ID'), header{end} = fn{i}; end
        C(:,end+1) = cell(size(C,1),1); %#ok<AGROW>
        C{1,end} = header{end};
        j = numel(header);
    end
    colOf(i) = j;
end
names = cellfun(@(x) char(string(x)),C(2:end,keyCol),'UniformOutput',false);
r = find(strcmp(names,rowKey),1);
if isempty(r), r = size(C,1)+1; C(r,:) = cell(1,size(C,2)); else, r = r+1; end
C(r,:) = cell(1,size(C,2));                           % clear the old values of this row
C{r,keyCol} = rowKey;
for i = 1:numel(fn)
    v = summary.(fn{i});
    if isstring(v), v = char(v); end
    if isnumeric(v) && isnan(v), v = []; end
    C{r,colOf(i)} = v;
end
writecell(C,xlsxPath,'Sheet',1);

% ---------- matrices sheet ----------
sn = regexprep(rowKey,'[\[\]\:\*\?\/\\]','_');
sn = sn(1:min(31,numel(sn)));
if any(strcmpi(sheetnames(xlsxPath),sn))
    old = readcell(xlsxPath,'Sheet',sn);
    writecell(cell(size(old)),xlsxPath,'Sheet',sn);   % blank out leftovers from a previous, larger export
end
writetable(dataTbl,xlsxPath,'Sheet',sn);
fprintf('  Analysis dataset: row "%s" updated in sheet 1; matrices written to sheet "%s"\n',rowKey,sn);
end

function x = clean(x)
if ismissing_(x), x = []; end
end

function tf = ismissing_(x)
tf = false;
if isa(x,'missing'), tf = true;
elseif isnumeric(x) && isscalar(x) && isnan(x), tf = true;
end
end
