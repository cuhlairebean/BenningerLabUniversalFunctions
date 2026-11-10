function t = getFrameTimes(R, howmanychannel)
%GETFRAMETIMES  Per-frame time stamps (seconds) from the Bio-Formats metadata of a .czi.
%   t = getFrameTimes(R, howmanychannel)
%   R is the output of bfopen. Uses the first channel's plane of each time point.
%   Returns [] if the file has no usable time stamps.
t = [];
try
    meta = R{1,4};
    nPlanes = meta.getPlaneCount(0);
    dt = nan(nPlanes,1);
    for p = 0:nPlanes-1
        v = meta.getPlaneDeltaT(0,p);
        if ~isempty(v), dt(p+1) = double(v.value()); end
    end
    dt = dt(1:howmanychannel:end);
    if any(~isnan(dt)), t = dt; end
catch
    t = [];
end
end
