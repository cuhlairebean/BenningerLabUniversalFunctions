function [DataOut]=Run_Nosilentcell(data, samefile, cafirst, threshhold)
%% LOADING IMAGES
DataOut = data;
if length(fieldnames(data))<3
    mask=[];
    bkgcell=[];   
else
    try
        mask=data.mask;
    catch
        mask=[];
    end
    
    try
        bkgcell=data.bkgcell;
    catch
        bkgcell=[];
    end
    
    try
        scthresh=data.scthresh;
    catch
        scthresh=[];
    end
    
    try
        bkgthresh=data.BKGthresh;
    catch
        bgthresh=[];
    end
end

R=bfopen(data.Location);
pics=R{1};
pics=pics(:,1);
pn = length(pics);

for i=1:pn
    IMG(:,:,i)=pics{i};
end

try
    for i=1:length(pics)
        T(i)=R{4}.getPlaneDeltaT(0, i-1).value;
    end
catch
    T=0:0.5:pn;
end
T = double(T);
clear pics

Images=double(IMG);

if samefile
<<<<<<< Updated upstream
    if cafirst
        Images = Images(:,:,1:2:end);
    else
        Images = Images(:,:,2:2:end);
    end
=======
    %Images = Images(:,:,whichchannelca,2:end);
    Image = Images(:,:,1:2:end);
>>>>>>> Stashed changes
end
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
imagesRaw=images;

h = fspecial('average', 5);
for i=1:sz
    images(:,:,i)=imfilter(imagesRaw(:,:,i),h,'replicate').*mask;
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
    dat = detrend(dat)/mean(Mean);%/mean(dat);
    
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

%%%% using the silent/background cell threshold to remove areas with NO signal %%%
SCMap=peakAmpMap.*bkgcell;
SCMap(~logical(SCMap))=nan;
scthresh = nanmean(nanmean(SCMap));

[NoSignalMask, ~]=RemovingAreaswNoSignal2(AverageIntensityMap,bkgthresh, peakAmpMap, scthresh);

SCMap=peakAmpMap.*NoSignalMask;
SCMap(~logical(SCMap))=nan;
scthresh = nanmedian(nanmedian(SCMap));

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

%LogicalMap(SmallAreaMap)=0;
Area=length(nonzeros(LogicalMap));
LogicalMap = logical(LogicalMap);

figure
overlay = imoverlay(mat2gray(RawImg), LogicalMap, 'colormap', 'jet', 'facealpha', .75, 'ZeroAlpha', 0);
imshow(overlay)

%% CALCULATING ACTIVITY and GENERATING ACTIVITY MAPS
ActMap=zeros(sx,sy);
%Determine the active area using the bachground threshhold
ActMap(peakAmpMap>=threshhold*scthresh)=1;
ActMapNoGFP=ActMap;

ActMap(~LogicalMap)=0;
ActMask=logical(ActMap);
ActMap(~ActMap)=NaN;

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

%% SAVING ACTIVITY VARIABLES
%set a breakpoint and look at graphs
DataOut.ActMap=ActMap;
DataOut.ActMapValues=CoorMap;
DataOut.CorrMap=coeffcorrMap;
DataOut.Location=data.Location;
DataOut.NoSignalMask=NoSignalMask;
DataOut.LogicalMap=LogicalMap;
DataOut.mask=mask;

%% Matt's correlation analysis
%[Dat]=HIanalyzeU(Images,T,0,LogicalMap);
CorrPeakAmp = zeros(sx,sy);
if (max(max(ActMask))>0)
    [Dat]=HIanalyzeU(imagesRaw,T,0,ActMask);
    Dat.Location=data.Location;
    Dat.Image=imagesRaw(:,:);
    Dat.ActMask = ActMask;
    Dat.T = T;
    Dat=CoordinatedAreas(Dat);
    
    CorrPeakAmp = Dat.Img();
    CorrPeakAmp = CorrPeakAmp;
    CorrPeakAmp(~LogicalMap)=0;
    stats = [];
    STATS=regionprops(CorrPeakAmp,'area');
    for j=1:length(STATS)
        stats(j)=STATS(j).Area;
    end
    if isempty(stats)
        maxcorrarea = 0;
        disp('maxcorr is zero')
    elseif length(max(stats))>1
        maxes = max(stats);
        maxcorrarea = maxes(1);
        disp('maxcorr - 2 areas with same size')
    else
        maxcorrarea = max(stats);
    end
else
    maxcorrarea = 0;
    Dat.AMCA = 0;
end
DataOut.CorrPeakAmpArea = maxcorrarea/Area;

figure()
imagesc(CorrPeakAmp);
colorbar;

disp('End of analysis.')
%sprintf('Active Area is %.4f %%', percentActiveArea)
end