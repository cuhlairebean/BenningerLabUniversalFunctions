function [DataOut]=GetCorrelationMapInterp(data,cachannel,howmanychannel,zstacks, starttime, endtime)
%% LOADING IMAGES
addpath('C:\Users\dwuletj\Documents\GitHub\UniversalCode\')
DataOut = data;

R=bfopen(DataOut.Location);
omeMeta = R{1, 4};
stackSizeX = omeMeta.getPixelsSizeX(0).getValue(); % image width, in pixels
stackSizeY = omeMeta.getPixelsSizeY(0).getValue(); % image height, in pixels

voxelSizeX = omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER); % in µm
voxelSizeXdouble = voxelSizeX.doubleValue(); % The numeric value represented by this object after conversion to type double
voxelSizeY = omeMeta.getPixelsPhysicalSizeY(0).value(ome.units.UNITS.MICROMETER); % in µm
voxelSizeYdouble = voxelSizeY.doubleValue(); % The numeric value represented by this object after conversion to type double

PixelSize = voxelSizeXdouble; % in um
clear R

%% full correlation map
if starttime == -1
    starttime = 1;
end
if endtime == -1
    endtime = size(DataOut.images,3);
end
T=starttime:endtime;
eachstack =  1.1542;
realTime = 1:eachstack:(endtime*zstacks*eachstack);
%T=starttime:eachstack*zstacks:(endtime*zstacks*eachstack);

newimages = [];
newimagesint = [];
newactmask = [];
newlogicalmap = [];
ActiveAreaStacks = 0;
AreaStacks = 0;

for zz=1:zstacks
    newimagesint = [newimagesint; DataOut.interpolatedimages(:,:,:,zz)];
    newimages = [newimages; DataOut.images(:,:,:,zz)];
    newactmask = [newactmask; DataOut.ActMask(:,:,zz)];
    newlogicalmap = [newlogicalmap; DataOut.LogicalMap(:,:,zz)];
    ActiveAreaStacks = ActiveAreaStacks + DataOut.ActiveArea(zz);
    AreaStacks = AreaStacks + DataOut.Area(zz);
end
[totalmaxcorrarea Dat] = getcoordarea(newactmask, newlogicalmap, newimages,T);
DataOut.CorrPeakAmpAreaStack = totalmaxcorrarea/AreaStacks;
DataOut.ActiveAreaStack = ActiveAreaStacks/AreaStacks;
corrpeakmap = Dat.TopCorrPeakAmp;

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

[~, TopcoordDat] = getcoordarea(mostcoordinated, newlogicalmap, newimagesint,realTime);
    
    phaseMap = TopcoordDat.phaseMap;
%     figure
%     imagesc(phaseMap);
    maxphase = max(max(phaseMap));
    minphase = min(min(phaseMap));
    difphase = maxphase-minphase;
    DataOut.difphase = difphase;
    %f = TopcoordDat.Freq
    avgfreq  = TopcoordDat.AvgFreq;
    period = TopcoordDat.period;
    
    %%need to do in same frame or across 3D
    %%calculate microns
    %then distance per time
    [rmax,cmax] = find(phaseMap==maxphase);
    [rmin,cmin] = find(phaseMap==minphase);
    
    currentdist = 0;
    x = [cmax rmax];
    y = [cmin rmin];
    for xx = 1:size(x,1)
        for yy =1:size(y,1)
            newx = x(xx,:);
            newy = y(yy,:);
            newdist = norm(newx-newy)*PixelSize;
            if newdist > currentdist
                DataOut.x = newx;
                DataOut.y = newy;
                DataOut.dist = newdist;
                currentdist =newdist;
            end
        end
    end

DataOut.avgdifphase = max(DataOut.difphase);
DataOut.allspeed = DataOut.dist./DataOut.difphase;
%DataOut.avgdist = DataOut.dist(find(max(DataOut.difphase)));
DataOut.maxspeed = max(DataOut.allspeed);
DataOut.maxdifphase = DataOut.difphase(find(DataOut.allspeed == DataOut.maxspeed));

disp('End of analysis.')
end