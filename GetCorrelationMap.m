function [DataOut]=GetCorrelationMap(data,cachannel,howmanychannel,zstacks, starttime, endtime)
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
newimages = [];
newimagesfilt = [];
newactmask = [];
newlogicalmap = [];
ActiveAreaStacks = 0;
AreaStacks = 0;

sx = size(DataOut.images,1);

for zz=1:zstacks
    newimages = [newimages; DataOut.images(:,:,:,zz)];
    %newimagesfilt = [newimagesfilt; DataOut.imagesfilt(:,:,:,zz)];
    newactmask = [newactmask; DataOut.ActMask(:,:,zz)];
    newlogicalmap = [newlogicalmap; DataOut.LogicalMap(:,:,zz)];
    ActiveAreaStacks = ActiveAreaStacks + DataOut.ActiveArea(zz);
    AreaStacks = AreaStacks + DataOut.Area(zz);
end
[totalmaxcorrarea Dat] = getcoordarea(newactmask, newlogicalmap, newimages,T);
DataOut.CorrPeakAmpAreaStack = totalmaxcorrarea/AreaStacks
DataOut.ActiveAreaStack = ActiveAreaStacks/AreaStacks;
DataOut.xcovmap = Dat.corrMap;
corrpeakmap = Dat.TopCorrPeakAmp;

% [totalmaxcorrareafilt Dat] = getcoordarea(newactmask, newlogicalmap, newimagesfilt,T);
% DataOut.CorrPeakAmpAreaStack = totalmaxcorrareafilt/AreaStacks;

close all

savefilename = DataOut.Location;
savefilename = regexprep(savefilename,'.xml','');


%% get phase for only coordinated areas
mostcoordinated = zeros(size(corrpeakmap));
%mostcoordinated(find(corrpeakmap==1))=1;

STATS=regionprops(corrpeakmap,'area');
for j=1:length(STATS)
    stats(j)=STATS(j).Area;
end
statsordered = sort(stats, 'descend');
r = find(stats==statsordered(1));
mostcoordinated(corrpeakmap==r)=1;

[~, TopcoordDat] = getcoordarea(mostcoordinated, newlogicalmap, newimages,T);

phaseMap = TopcoordDat.phaseMap;
figure
imagesc(phaseMap);
maxphase = max(max(phaseMap));
minphase = min(min(phaseMap));
difphase = maxphase-minphase;
%f = TopcoordDat.Freq
avgfreq  = TopcoordDat.AvgFreq;
period = TopcoordDat.period;

%%need to do in same frame or across 3D
%%calculate microns
%then distance per time
[rmax,cmax] = find(phaseMap==maxphase)
[rmin,cmin] = find(phaseMap==minphase)

x = [cmax rmax];
y = [cmin rmin];
norm(x-y)

corrpeakmap(~logical(corrpeakmap))=NaN;
st=1;
%% SET MAP COLORS
for zz=1:zstacks 
    RawImg = sum(DataOut.images(:,:,:,zz),3);
    pause(2)
    %figure()
    imoverlayNurin(mat2gray(RawImg),corrpeakmap(st:st+sx-1,:),[],[],'parula',0.3);
    caxis([0 5]);
    %colorbar
    %title('Correlated Areas')
    imagewd = getframe(gcf);
    imwrite(imagewd.cdata, [savefilename num2str(zz) '.tiff'], 'Resolution',600);
    st=st+sx;
end

disp('End of analysis.')
end