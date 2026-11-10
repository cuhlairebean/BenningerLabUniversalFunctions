function ratio = campariCellRatios(CellMask, redImg, greenImg, bg)
%CAMPARICELLRATIOS  Background-subtracted red:green ratio of every cell (same formula as selectCampariThreshold).
%
%   ratio = campariCellRatios(CellMask, redImg, greenImg, bg)
%
%   CellMask  label matrix (cell j = j) on the time-series image grid (e.g. 512x512)
%   redImg, greenImg  channels of the high-res _pre.czi (any size; the masks are stretched onto them)
%   bg        [greenBackground redBackground]
%   ratio     1 x max(label) vector: (mean red in cell - redBg) / (mean green in cell - greenBg), NaN if the cell is empty
n = max(CellMask(:));
ratio = nan(1,n);
if n == 0, return; end
redImg = double(redImg); greenImg = double(greenImg);
M = imresize(CellMask,[size(redImg,1) size(redImg,2)],'nearest');
for j = 1:n
    m = (M == j);
    if ~any(m(:)), continue; end
    ratio(j) = (mean(redImg(m)) - bg(2)) / (mean(greenImg(m)) - bg(1));
end
end
