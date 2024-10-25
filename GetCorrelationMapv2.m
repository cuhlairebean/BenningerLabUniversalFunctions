function [DataOut]=GetCorrelationMapv2(data,zstacks,PixelSize)
%% LOADING IMAGES
addpath('C:\Users\dwuletj\Documents\GitHub\UniversalCode\')
DataOut = data;

%% full correlation map
% if starttime == -1
%     starttime = 1;
% end
% if endtime == -1
%     endtime = size(DataOut.images,3);
% end
% T=starttime:endtime;
T=1:size(DataOut.images,3);

newimages = [];
newimagesfilt = [];
newactmask = [];
newlogicalmap = [];
ActiveAreaStacks = 0;
AreaStacks = 0;
newphaseofcorrarea = [];

sx = size(DataOut.images,1);

for zz=1:zstacks
    newimages = [newimages; DataOut.images(:,:,:,zz)];
    %newimagesfilt = [newimagesfilt; DataOut.imagesfilt(:,:,:,zz)];
    newactmask = [newactmask; DataOut.ActMask(:,:,zz)];
    newlogicalmap = [newlogicalmap; DataOut.LogicalMap(:,:,zz)];
    ActiveAreaStacks = ActiveAreaStacks + DataOut.ActiveArea(zz);
    AreaStacks = AreaStacks + DataOut.Area(zz);
    newphaseofcorrarea = [newphaseofcorrarea; DataOut.PhaseMap3D(:,:,zz)];
end
DataOut.ActiveAreaStack = ActiveAreaStacks/AreaStacks;

[totalmaxcorrarea, Dat] = getcoordarea(newactmask, newlogicalmap, newimages,T);
DataOut.CorrPeakAmpAreaStack = totalmaxcorrarea/AreaStacks;
corrpeakmap = Dat.TopCorrPeakAmp;
DataOut.topcorrpeakmap = corrpeakmap;
DataOut.allcorrpeakmap = Dat.CorrPeakAmp;
DataOut.maxcorrarea_num = Dat.maxcorrarea_num;

largestarea = corrpeakmap;
largestarea(largestarea ~= Dat.maxcorrarea_num)=0;

newphaseofcorrarea = largestarea.*newphaseofcorrarea;
figure;
imagesc(newphaseofcorrarea);
title('orig')
colorbar




% disp('select edge')
% bf=imfreehand();
% useMask=createMask(bf);
% onepoint=useMask;
% stats1 = regionprops('table',onepoint,'Centroid')
% centers1 = stats1.Centroid;
% 
% disp('select edge')
% bf=imfreehand();
% useMask=createMask(bf);
% otherpoint=useMask;
% stats2 = regionprops('table',otherpoint,'Centroid')
% centers2 = stats2.Centroid;
% DataOut.distofcorrarea = norm(centers1-centers2)*PixelSize

% newphaseofcorrarea = smoothimage(newphaseofcorrarea, 5);
% figure;
% imagesc(newphaseofcorrarea);
% title('smooth')
% colorbar

% fun = @(block_struct) mean(block_struct.data(:));
% newphaseofcorrareabloc = blockproc(newphaseofcorrarea, [10 10], fun);
% 
% newphaseofcorrareabloc(newphaseofcorrareabloc==0)=NaN;
% 
% figure;
% imagesc(newphaseofcorrareabloc);
% title('phase blocked');
% colorbar

% maxphase = max(newphaseofcorrarea,[],'all');
% minphase = min(newphaseofcorrarea,[],'all');

%close all

% savefilename = DataOut.Location;
% savefilename = regexprep(savefilename,'.xml','');
% 
% DataOut.difphase = maxphase-minphase;
% 
% [rmax,cmax] = find(newphaseofcorrareabloc==maxphase);
% [rmin,cmin] = find(newphaseofcorrareabloc==minphase);
% 
% x = [cmax rmax];
% y = [cmin rmin];
% DataOut.distminmax = norm(x-y)*PixelSize;

% DataOut.allspeed = DataOut.dist./DataOut.difphase
%DataOut.allspeed = freq*DataOut.dist./DataOut.difphase
% DataOut.maxdist = DataOut.dist(find(DataOut.difphase == max(DataOut.difphase)));
% DataOut.speed = DataOut.allspeed(find(DataOut.difphase == DataOut.maxdifphase));


% corrpeakmap(~logical(corrpeakmap))=NaN;
% st=1;
% %% SET MAP COLORS
% for zz=1:zstacks
%     RawImg = sum(DataOut.images(:,:,:,zz),3);
%     pause(2)
%     %figure()
%     imoverlayNurin(mat2gray(RawImg),corrpeakmap(st:st+sx-1,:),[],[],'parula',0.3);
%     caxis([0 5]);
%     %colorbar
%     %title('Correlated Areas')
%     imagewd = getframe(gcf);
%     imwrite(imagewd.cdata, [savefilename num2str(zz) '.tiff'], 'Resolution',600);
%     st=st+sx;
% end

disp('End of analysis.')
end