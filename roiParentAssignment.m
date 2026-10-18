function [parentIdx, overlapFrac] = roiParentAssignment(roiLabel, parentLabel)
%ROIPARENTASSIGNMENT  Which parent region (e.g. islet #) does each ROI fall in?
%   roiLabel    label matrix of ROIs (1..N)
%   parentLabel label matrix of parents (1..M), same size
%   parentIdx   Nx1, parent number holding the majority of each ROI's pixels
%               (0 = ROI not inside any parent)
%   overlapFrac Nx1, fraction of the ROI's pixels inside that parent
N = max(roiLabel(:));
parentIdx = zeros(N,1); overlapFrac = zeros(N,1);
for n = 1:N
    m = roiLabel == n;
    tot = nnz(m);
    pix = parentLabel(m); pix = pix(pix>0);
    if tot == 0 || isempty(pix), continue; end
    parentIdx(n) = mode(pix);
    overlapFrac(n) = nnz(pix==parentIdx(n))/tot;
end
end
