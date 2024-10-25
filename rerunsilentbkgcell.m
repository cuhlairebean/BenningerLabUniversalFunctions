% bkgcell = []
% silentcell = []

%% APPLYING "SILENT CELL" THRESHOLD

%%%%%%Selecting "background cell" from peak amp map%%%%%%%%%%
AverageIntensityMap = mean(images, 3).*mask;
figure; imagesc(AverageIntensityMap);
colorbar;
if isempty(bkgcell)
    figure;
    %imagesc(AverageIntensityMap/(mean(mean(AverageIntensityMap)))+peakAmpMap/(mean(mean(peakAmpMap))))
    imagesc(AverageIntensityMap)
    disp('select background cell')
    bkgcell=createMask(imfreehand);
else
    figure
    overlay = imoverlay(mat2gray(RawImg), bkgcell, 'colormap', 'jet', 'facealpha', .75, 'ZeroAlpha', 0);
    imshow(overlay)
end
BKGMap=AverageIntensityMap.*bkgcell;
BKGMap(~logical(BKGMap))=nan;
bkgthresh=nanmax(nanmax(BKGMap));

DataOut.BKGMap=BKGMap;
DataOut.BKGthresh=bkgthresh;
DataOut.bkgcell=bkgcell;

%%%%%%%%%%%%silentcell%%%%%%%%%%%%%%%%%
silentcell = bkgcell;
if isempty(silentcell)
    figure;
    imagesc(peakAmpMap)
    disp('select silent cell')
    silentcell=createMask(imfreehand);
else
    figure
    overlay = imoverlay(mat2gray(RawImg), silentcell, 'colormap', 'jet', 'facealpha', .75, 'ZeroAlpha', 0);
    imshow(overlay)
end

%%%% using the silent/background cell threshold to remove areas with NO signal %%%
SCMap=peakAmpMap.*bkgcell;
SCMap(~logical(SCMap))=nan;
scthresh = nanmean(nanmean(SCMap));

%[NoSignalMask, NoSignalArea]=RemovingAreaswNoSignal(peakAmpMap,bgthresh);
%[NoSignalMask2, NoSignalArea2]=RemovingAreaswNoSignal(AverageIntensityMap,bgthresh2);
%NoSignalMask = NoSignalMask2;
[NoSignalMask, ~]=RemovingAreaswNoSignal2(AverageIntensityMap,bkgthresh, peakAmpMap, scthresh);

SCMap=peakAmpMap.*NoSignalMask;
SCMap(~logical(SCMap))=nan;
scthresh = nanmax(nanmax(SCMap));

DataOut.silentcell=silentcell;
DataOut.scthresh=scthresh;

figure
overlay = imoverlay(mat2gray(RawImg), NoSignalMask, 'colormap', 'jet', 'facealpha', .75, 'ZeroAlpha', 0);
imshow(overlay)
%% %%%%%% Final Mask to use%%%%%%%%%%%%%%%%%%%%%%%%
LogicalMap=zeros(sx,sy);
LogicalMap(logical(mask))=1;
LogicalMap(NoSignalMask)=0;
TotalArea=length(nonzeros(LogicalMap));
DataOut.TotalArea = TotalArea;

if useGFPasMaskBoolean
    %notdtom area
    LogicalMapNoGFP = LogicalMap;
    LogicalMapNoGFP(GFPlogical)=0;
    LogicalMapNoGFP = logical(LogicalMapNoGFP);
    NoGFPArea = length(nonzeros(LogicalMapNoGFP));
    %tdtom area
    LogicalMap(~GFPlogical)=0;
    figure
    overlay = imoverlay(mat2gray(GFPimageorig), LogicalMap, 'colormap', 'jet', 'facealpha', .75, 'ZeroAlpha', 0);
    imshow(overlay)
end

%LogicalMap(SmallAreaMap)=0;
Area=length(nonzeros(LogicalMap));
LogicalMap = logical(LogicalMap);

figure
overlay = imoverlay(mat2gray(RawImg), LogicalMap, 'colormap', 'jet', 'facealpha', .75, 'ZeroAlpha', 0);
imshow(overlay)


if useGFPasMaskBoolean
    %GFPArea = size(nonzeros(GFPlogical),1);
    %DataOut.GFPArea = GFPArea;
    DataOut.GFPPercent = Area/TotalArea;
    DataOut.NotGFPPercent = length(nonzeros(LogicalMapNoGFP))/TotalArea;
    %Area = Area - GFPArea;
end


%% CALCULATING ACTIVITY and GENERATING ACTIVITY MAPS
ActMap=zeros(sx,sy);
%Determine the active area using the bachground threshhold

ActMap(peakAmpMap>=3)=1;
ActMapNoGFP=ActMap;

ActMap(~LogicalMap)=0;
ActMask=logical(ActMap);
ActMap(~ActMap)=NaN;

if useGFPasMaskBoolean
    ActMapNoGFP(~LogicalMapNoGFP)=0;
    ActMaskNoGFP = logical(ActMapNoGFP);
    ActMapNoGFP(~ActMapNoGFP)=NaN;
end


%ActMap(peakAmpMap<2.5*mean([bgthresh, scthresh]))=0;
%%%%%% Remove small areas from the activity map%%%%%%%%%%%%%%%%%%%%%%%%%
%%Reason: Small areas are probably not cells
% RemoveSmallArea = 0;
% if (RemoveSmallArea)
%     SmallAreaMap = RemovingSmallAreas(ActMap, 100);
%     ActMap(logical(SmallAreaMap))=0;
%     SmallAreaMap=logical(SmallAreaMap);
%
%     pause(2)
%     figure;
%     imagesc(SmallAreaMap);
%     colormap('jet')
%     title('Small Area Map for Activity');
% end

%% SET MAP COLORS
colors = [
    1 0 1 % First element = purple
    0 0 1 % blue
    0 1 0 % green
    1 1 0 % yellow
    1 .65 0 % orange
    1 0 0]; % red

n = 256; % size of new color map
m = size(colors,1);
t0 = linspace(0,1,m)';
tt = linspace(0,1,n)';
r = interp1(t0,colors(:,1),tt);
g = interp1(t0,colors(:,2),tt);
b = interp1(t0,colors(:,3),tt);
cmap2 = [r,g,b];


%%%%%%Use the CoorMap to overlay the activity map to put a number to activityremoving regions in the no signal mask, outside of the activity mask, and the image mask
CoorMap=coeffcorrMap;
CoorMap(~ActMask)=0;

%CoorMap(CoorMap<.15)=0;  %%Another threshhold that is useful sometimes especially with movement
ActiveArea=size(nonzeros(CoorMap),1);
DataOut.RatioActive=ActiveArea/(Area);%% CALCULATING PERCENT ACTIVE AREA
DataOut.ActiveArea = ActiveArea;
DataOut.Area = Area;
disp(['Ratio Active ' num2str(DataOut.RatioActive)]) %print active are

%% PLOTTING ACTIVITY MAP
CoorMap(~logical(CoorMap))=NaN;
pause(2)
imoverlayNurin(mat2gray(imagesRaw(:,:,1)),CoorMap,[0, 1],[],cmap2,0.3);
colorbar
caxis([0 1])
title('Activity Map')

if useGFPasMaskBoolean
    CoorMapNoGFP=coeffcorrMap;
    CoorMapNoGFP(~ActMaskNoGFP)=0;
    %CoorMap(CoorMap<.15)=0;  %%Another threshhold that is useful sometimes
    ActiveAreaNoGFP=size(nonzeros(CoorMapNoGFP),1);
    DataOut.RatioActiveNoGFP=ActiveAreaNoGFP/(NoGFPArea);%% CALCULATING PERCENT ACTIVE AREA
    CoorMapNoGFP(~logical(CoorMapNoGFP))=NaN;
    pause(2)
    imoverlayNurin(mat2gray(imagesRaw(:,:,1)),CoorMapNoGFP,[0, 1],[],cmap2,0.3);
    colorbar
    caxis([0 1])
    title('Activity Map No Td Tom')
    DataOut.ActMapNoGFP=ActMapNoGFP;
    DataOut.ActMapValuesNoGFP=CoorMapNoGFP;
end

%% SAVING ACTIVITY VARIABLES
%set a breakpoint and look at graphs
DataOut.ActMap=ActMap;
DataOut.ActMapValues=CoorMap;
DataOut.CorrMap=coeffcorrMap;
DataOut.Location=data.Location;
DataOut.NoSignalMask=NoSignalMask;
DataOut.LogicalMap=LogicalMap;
DataOut.mask=mask;