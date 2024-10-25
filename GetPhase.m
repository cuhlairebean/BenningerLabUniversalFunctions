function [DataOut]=GetPhase(data,cachannel,howmanychannel,zstacks, starttime, endtime)
%% LOADING IMAGES
addpath('C:\Users\dwuletj\Documents\GitHub\UniversalCode\')
DataOut = data;

%% full correlation map
if starttime == -1
    starttime = 1;
end
if endtime == -1
    endtime = size(DataOut.images,3);
end
T=starttime:endtime;

for zz=1:zstacks
    newimages = DataOut.images(:,:,:,zz);
    newactmask = DataOut.ActMask(:,:,zz);
    newlogicalmap = DataOut.LogicalMap(:,:,zz);
    ActiveAreaStacks = DataOut.ActiveArea(zz);
    AreaStacks = DataOut.Area(zz);
    [~, Dat] = getcoordarea(newactmask, newlogicalmap, newimages,T);
    phaseMap = Dat.phaseMap;
    maxphase = max(max(phaseMap));
    minphase = min(min(phaseMap));
    difphase(zz,1) = maxphase-minphase;
end
DataOut.difphase = difphase;
DataOut.maxphase = max(difphase)

% %% get phase for only coordinated areas
% mostcoordinated = zeros(size(corrpeakmap));
% %mostcoordinated(find(corrpeakmap==1))=1;
%
% STATS=regionprops(corrpeakmap,'area');
% for j=1:length(STATS)
%     stats(j)=STATS(j).Area;
% end
% statsordered = sort(stats, 'descend');
% r = find(stats==statsordered(1));
% mostcoordinated(corrpeakmap==r)=1;
%
% [~, TopcoordDat] = getcoordarea(mostcoordinated, newlogicalmap, newimages,T);
%
% phaseMap = TopcoordDat.phaseMap;
% maxphase = max(max(phaseMap))
% minphase = min(min(phaseMap))
% difphase = maxphase-minphase
% %f = TopcoordDat.Freq
% avgfreq  = TopcoordDat.AvgFreq
% period = TopcoordDat.period


disp('End of analysis.')
end