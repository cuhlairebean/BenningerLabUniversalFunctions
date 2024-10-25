close all 
%mask=[];
%bkgcell=[];
bkgcell=[];
 
%% CODE TO MODIFY TCs
%Can adjust the time frame of images here%
st=1;
ed=size(Images,3);
T=st:ed;
%t=st:ed; %for duty cycle

images=Images(:,:,st:ed);
DataOut.images=images;
%Images=double(images);
RawImg=images(:,:,1);
MaxIntensityMap = max(images, [],3);

%% DECLARING IMAGE PROPERTIES

sx=size(images,1);
sy=size(images,2);
sz=size(images,3);

%% MASKING IMAGES

if isempty(mask)
    figure
    imshow(mat2gray(MaxIntensityMap));
    numAreas=1;
    mask=zeros(sx,sy);
    % level=graythresh(mat2gray(RawImg));
    % mask=im2bw(mat2gray(RawImg),level);
    disp('select islet')
    for i=1:numAreas
        bf=imfreehand();
        useMask=createMask(bf);
        mask=mask+useMask;
    end
else
    figure
    overlay = imoverlay(mat2gray(RawImg), mask, 'colormap', 'jet', 'facealpha', .75, 'ZeroAlpha', 0);
    imshow(overlay);
end
DataOut.mask=mask;
%labelbw=bwlabel(mask,4);
%numAreas=max(labelbw(:));
%currArea=1;
imagesRaw=images;

h = fspecial('average', 5);
for i=1:sz
    images(:,:,i)=imfilter(imagesRaw(:,:,i),h,'replicate').*mask;
    %images(:,:,i)=images(:,:,i).*mask; %Nurin
    %images(:,:,i)=medfilt2(imagesRaw(:,:,i),[5 5]).*mask;
end

%% CALCULATING MAPS
timeVec=T;
Fs=mean(diff(timeVec));
Fs=1/Fs;

Mean = squeeze(mean(mean(images,1),2));

%oscMap=zeros(sx,sy);
peakMap=zeros(sx,sy);
troughMap=zeros(sx,sy);
timeMat=zeros(sx,sy);
peakAmpMap=zeros(sx,sy);
corrMap=nan(sx,sy);
coeffcorrMap=zeros(sx,sy);
savedavgvalues=zeros(sx,sy,sz);
stddevMap = zeros(sx,sy);
t=1:sz;

[x,y]=find(mask==1);

for i=1:length(x)
    
    xst = x(i)-5;
    xed = x(i)+5;
    yst = y(i)-5;
    yed = y(i)+5;
    if xst < 1
        xst = 1;
    end
    if xed > 512
        xed = 512;
    end
    if yst < 1
        yst = 1;
    end
    if yed > 512
        yed = 512;
    end
    
    dat = images(xst:xed,yst:yed,:);
    % dat = images(x(i)-3:x(i)+3,y(i)-3:y(i)+3,:); %Nurin
    
    dat = mean(dat,2);
    dat = mean(dat,1);
    dat = dat(:);
    
    f = polyfit(t',dat,2);
    feval = polyval(f,t');
    dat2 = dat./feval;
    %dat = dat-feval;
    dat = detrend(dat)/mean(Mean);%/mean(dat);
%     figure
%     plot(dat)
    
    savedavgvalues(x(i), y(i), :) = dat2;
    
    %%try standard deviation
    stddevMap(x(i),y(i))=std(dat2);
    [~,peakAmp1]=peakfinder(dat2,(max(dat2)-min(dat2))/5);
    %dat=detrend(dat);
    
    [peakLoc,peakAmp]=peakfinder(dat,(max(dat)-min(dat))/5);
    [trough,troughAmp]=peakfinder(dat,(max(dat)-min(dat))/5,[],-1);
    numpeaks=size(peakLoc,1);
    peakMap(x(i),y(i))=numpeaks;
    troughMap(x(i),y(i))=size(trough,1);
    if length(peakLoc)>1
        timeMat(x(i),y(i))=mean(diff(peakLoc,1));
        peakAmpMap(x(i),y(i))=(mean(peakAmp)-mean(troughAmp));
    end
    if length(peakLoc)==1
        timeMat(x(i),y(i))=peakLoc;
        peakAmpMap(x(i),y(i))=peakAmp-mean(troughAmp);
    end
    if isempty(peakLoc)
        timeMat(x(i),y(i))=0;
        peakAmpMap(x(i),y(i))=0;
    end
    
    [c, lags]=xcov(dat,Mean,'coeff');
    [maxCV, maxCL]=max(c);
    
    if ~isnan(dat(1))
        corrMap(x(i),y(i))=lags(maxCL)/Fs; %phase
        coeffcorrMap(x(i),y(i))= maxCV; %coefficient
    end
    
end
DataOut.xInfo=length(x);
DataOut.PeakAmpMap=peakAmpMap;
DataOut.StdDevMap = stddevMap;

figure;
imagesc(peakAmpMap);
colorbar;
title('Peak Amp Map');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%ACTIVITY CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% APPLYING "SILENT CELL" THRESHOLD

%%%%%%Selecting "background cell" from peak amp map%%%%%%%%%%
AverageIntensityMap = mean(images, 3).*mask;
figure; imagesc(AverageIntensityMap);
colorbar;
if isempty(bkgcell)
    figure;
    %imagesc(AverageIntensityMap/(mean(mean(AverageIntensityMap)))+peakAmpMap/(mean(mean(peakAmpMap))))
    imagesc(peakAmpMap)
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
scthresh = nanmedian(nanmedian(SCMap));

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

ActMap(peakAmpMap>=2*scthresh)=1; %also tried 2.3
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