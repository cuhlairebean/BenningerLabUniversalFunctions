function [DataOut]=GetCorrelationMapInterpolation(data,cachannel,howmanychannel,zstacks, starttime, endtime)
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

%% full correlation map
if starttime == -1
    starttime = 1;
end
if endtime == -1
    endtime = size(DataOut.images,3);
end
T=starttime:endtime;
%eachstack =  1.1542;
%realTime = 1:eachstack:(endtime*zstacks*eachstack);
%T=starttime:eachstack*zstacks:(endtime*zstacks*eachstack);

newimages = [];
newimagesint = [];
newimagesfilt = [];
newactmask = [];
newlogicalmap = [];
ActiveAreaStacks = 0;
AreaStacks = 0;

for zz=1:zstacks
    % newimagesint = [newimagesint; DataOut.interpolatedimages(:,:,:,zz)];
    newimages = [newimages; DataOut.images(:,:,:,zz)];
    newimagesfilt = [newimagesfilt; DataOut.imagesfilt(:,:,:,zz)];
    newactmask = [newactmask; DataOut.ActMask(:,:,zz)];
    newlogicalmap = [newlogicalmap; DataOut.LogicalMap(:,:,zz)];
    ActiveAreaStacks = ActiveAreaStacks + DataOut.ActiveArea(zz);
    AreaStacks = AreaStacks + DataOut.Area(zz);
end
[totalmaxcorrarea Dat] = getcoordarea(newactmask, newlogicalmap, newimages,T);
DataOut.CorrPeakAmpAreaStack = totalmaxcorrarea/AreaStacks;
DataOut.ActiveAreaStack = ActiveAreaStacks/AreaStacks;
corrpeakmap = Dat.TopCorrPeakAmp;

phaseMap1 = Dat.phaseMap;
phaseMap1(find(corrpeakmap==0))=0;

savefilename = DataOut.Location;
savefilename = regexprep(savefilename,'.xml','');
ztime = 9.371197817;

begin = 1;
for zz=1:zstacks
    phaseMap = phaseMap1(begin:begin+511,:);
    figure
    imagesc(phaseMap);
    colorbar
    
    mask=zeros(size(phaseMap));
    fullmask=zeros(size(phaseMap));
    
    disp('select islet')
    bf=imfreehand();
    useMask=createMask(bf);
    maxmask=mask+useMask;
    fullmask = mask+useMask;
    
    mask=zeros(size(phaseMap));
    bf=imfreehand();
    useMask=createMask(bf);
    minmask=mask+useMask;
    fullmask = fullmask +useMask.*2;
    
    figure
    imagesc(fullmask);
    
    phaseMapmax= phaseMap;
    phaseMapmin= phaseMap;
    
    phaseMapmax(find(maxmask==0))=nan;
    phaseMapmin(find(minmask==0))=nan;
    DataOut.phaseMapmax = phaseMapmax;
      DataOut.phaseMapmin = phaseMapmin;
    
    figure
    imagesc(phaseMapmax);
        
    figure
    imagesc(phaseMapmin);
    
    maxphase = nanmax(nanmax(phaseMapmax));
    minphase = nanmin(nanmin(phaseMapmin));
    minphase = nanmin(nanmin(phaseMapmax));
    difphase = maxphase-minphase;
    DataOut.difphase(zz) = difphase*ztime;
    
    STATS=regionprops(bwlabeln(imbinarize(fullmask)),'centroid');
    
    [rmax,cmax] = find(phaseMapmax==maxphase);
    [rmin,cmin] = find(phaseMapmax==minphase);
    
    x = [cmax rmax];
    y = [cmin rmin];
    
    newx = x(1,:);
    newy = y(1,:);
    size(x,1)
    size(y,1)
    currentdist = norm(newx-newy)/PixelSize;
     DataOut.x(zz,:) = newx;
    DataOut.y(zz,:) = newy;
    DataOut.dist(zz) = currentdist;
    totaldist = currentdist;
    count = 1;
    for xx = 1:size(x,1)
        for yy =1:size(y,1)
            newx = x(xx,:);
            newy = y(yy,:);
            newdist = norm(newx-newy)/PixelSize;
            totaldist = totaldist+newdist;
            count=count+1;
            if newdist < currentdist
                DataOut.x(zz,:) = newx;
                DataOut.y(zz,:) = newy;
                DataOut.dist(zz) = newdist;
                currentdist =newdist;
                totaldist = totaldist+newdist;
            end
        end
    end
    
%     newx = STATS(1).Centroid;
%     newy = STATS(2).Centroid;
%     newdist = norm(newx-newy)/PixelSize;
%     DataOut.x(zz,:) = newx;
%     DataOut.y(zz,:) = newy;
%     DataOut.dist(zz) = newdist;
    begin=begin+512;
end

DataOut.maxdifphase = max(DataOut.difphase);
DataOut.allspeed = DataOut.dist./DataOut.difphase
DataOut.maxdist = DataOut.dist(find(DataOut.difphase == max(DataOut.difphase)));
DataOut.speed = DataOut.allspeed(find(DataOut.difphase == DataOut.maxdifphase));

disp('End of analysis.')
end