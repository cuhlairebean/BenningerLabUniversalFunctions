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
% 
% mostcoordinated = zeros(size(corrpeakmap));
% %mostcoordinated(find(corrpeakmap==1))=1;
% STATS=regionprops(bwlabeln(imbinarize(round(phaseMap1,2))));
% STATS=regionprops(phaseMap1,'area');
% for j=1:length(STATS)
%     stats(j)=STATS(j).Area;
% end
% statsordered = sort(stats, 'descend');
% r = find(stats==statsordered(1));
% %mostcoordinated(find(corrpeakmap==r))=1;
% 
% 
% smallareamap = RemovingSmallAreas( phaseMap1, 200 );
% phaseMap1(find(smallareamap==1))=0;
% %mostcoordinated(corrpeakmap==r)=1;
% figure
% imagesc(phaseMap1);

% figure
% imagesc(corrpeakmap);

% [totalmaxcorrareafilt Dat] = getcoordarea(newactmask, newlogicalmap, newimagesfilt,T);
% DataOut.CorrPeakAmpAreaStack = totalmaxcorrareafilt/AreaStacks;
st = 1;
%close all

savefilename = DataOut.Location;
savefilename = regexprep(savefilename,'.xml','');

begin = 1;
for zz=1:zstacks
    % [~, TopcoordDat] = getcoordarea(mostcoordinated(begin:begin+511,:),DataOut.LogicalMap(:,:,zz), DataOut.images(:,:,:,zz),T);
    phaseMap = phaseMap1(begin:begin+511,:);
    maxphase = max(max(phaseMap));
    minphase = min(min(phaseMap));
    difphase = maxphase-minphase;
    DataOut.difphase(zz) = difphase;
    %f = TopcoordDat.Freq
    %     avgfreq  = TopcoordDat.AvgFreq;
    %     period = TopcoordDat.period;
    
    %%need to do in same frame or across 3D
    %%calculate microns
    %then distance per time
    [rmax,cmax] = find(phaseMap==maxphase);
    [rmin,cmin] = find(phaseMap==minphase);
    
    x = [cmax rmax];
    y = [cmin rmin];
    
    newx = x(1,:);
    newy = y(1,:);
    size(x,1)
    size(y,1)
    currentdist = norm(newx-newy)/PixelSize;;
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
    DataOut.totaldist(zz,:) = totaldist/count;
    %     DataOut.x(zz,:) = x;
    %     DataOut.y(zz,:) = y;
    %     DataOut.dist(zz) = norm(x-y)*PixelSize;
    begin=begin+512;
end
DataOut.maxdifphase = max(DataOut.difphase);
DataOut.allspeed = DataOut.dist./DataOut.difphase;
DataOut.maxdist = DataOut.dist(find(DataOut.difphase == max(DataOut.difphase)));
DataOut.meanallspeed = DataOut.totaldist./DataOut.difphase;
DataOut.meandist = DataOut.totaldist(find(DataOut.difphase == max(DataOut.difphase)));
count = 2;
%     while DataOut.maxdist < 2*PixelSize
%         phasesort = sort(DataOut.difphase);
%         DataOut.maxdifphase =phasesort(count);
%         DataOut.maxdist = DataOut.dist(find(DataOut.difphase == DataOut.maxdifphase));
%     end
DataOut.speed = DataOut.allspeed(find(DataOut.difphase == DataOut.maxdifphase))
DataOut.meanspeed = DataOut.meanallspeed(find(DataOut.difphase == DataOut.maxdifphase))

%% SET MAP COLORS
% for zz=1:zstacks
%     RawImg = sum(DataOut.images(:,:,:,zz),3);
%     pause(2)
%     %figure()
%     imoverlayNurin(mat2gray(RawImg),corrpeakmap(st:st+511,:),[],[],'parula',0.3);
%     caxis([0 5]);
%     %colorbar
%     %title('Correlated Areas')
%     imagewd = getframe(gcf);
%     imwrite(imagewd.cdata, [savefilename num2str(zz) '.tiff'], 'Resolution',600);
%     st=st+512;
% end
disp('End of analysis.')
end