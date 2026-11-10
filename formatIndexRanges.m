function txt = formatIndexRanges(ids)
% FORMATINDEXRANGES  Compress integers into ranges: [1 2 3 4 7 9 10] -> '1-4, 7, 9-10'
ids = unique(ids(:)');
if isempty(ids), txt = ''; return; end
br = [0 find(diff(ids) > 1) numel(ids)];
parts = cell(1,numel(br)-1);
for i = 1:numel(br)-1
    a = ids(br(i)+1); b = ids(br(i+1));
    if a == b, parts{i} = sprintf('%d',a); else, parts{i} = sprintf('%d-%d',a,b); end
end
txt = strjoin(parts,', ');
end
