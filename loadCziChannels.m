function stack = loadCziChannels(file)
%LOADCZICHANNELS  Read a .czi/.lsm/.tif with Bio-Formats into an HxWxC stack.
%   Each channel is the max projection across any Z / time planes.
%   Requires bfopen (Bio-Formats) on the path.
R = bfopen(file);
planes = R{1,1}(:,1);
nPlanes = numel(planes);
try
    nC = R{1,4}.getPixelsSizeC(0).getValue();
catch
    nC = 1;
end
nC = max(1,min(double(nC),nPlanes));
[h,w] = size(planes{1});
stack = zeros(h,w,nC,'single');
for c = 1:nC
    idx = c:nC:nPlanes;
    ch = single(planes{idx(1)});
    for k = 2:numel(idx)
        ch = max(ch,single(planes{idx(k)}));
    end
    stack(:,:,c) = ch;
end
end
