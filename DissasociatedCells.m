function [DataOut]=DissasociatedCells(data, GFP)
%% LOADING IMAGES
DataOut = data;
activitycutoff = 0.25;
allmasks = data.allmasks;
maxPeakThresh = 3;

try
    GFPbgMask = data.GFPbgMask;
catch
    GFPbgMask = [];
end

try
    CabgMask = data.CabgMask;
catch
    CabgMask = [];
end

%% %%%%%%%%%Open Calcium images%%%%%%%%%%%%%%%%%%%
R=bfopen(data.Location);
pics=R{1};
pics=pics(:,1);

% try R{4}.getPlaneDeltaT(0,2)
%
%     if ~isempty(R{4}.getPlaneDeltaT(0,0))
%         for i=1:length(pics)
%             T(i)=R{4}.getPl/aneDeltaT(0,i-1);
%         end
%     else
%         T=1:size(pics,3);
%     end
%
% catch
%     T=0:.5:length(pics);
% end
% T=double(T);

for i=1:length(pics)
    IMG(:,:,i)=pics{i};
end
clear pics;

Images=double(IMG(:,:,2:end));
RawImg=Images(:,:,2);
clear IMG;

%% %%%%%%%%%%%%%%%%%%%% Open GFP image %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R=bfopen(GFP.Location);
pics=R{1};
pics=pics(:,1);

for i=1:length(pics)
    IMG(:,:,i)=pics{i};
end
clear pics;

try
    GFPimages=double(IMG(:,:,2));
catch
    GFPimages=double(IMG(:,:,1));
end
GFPimagesRaw = GFPimages;

%% CODE TO MODIFY TCs
st=1;
ed=size(Images,3);
T=st:.5:ed;
images=Images(:,:,st:ed);
imagesRaw=images;
DataOut.images=images;

% % % %
%% DECLARING IMAGE PROPERTIES

sx=size(images,1);
sy=size(images,2);
sz=size(images,3);

%% Select Background GFP Image
allactivity=zeros(sx,sy);
labeledmask = zeros(sx,sy);

if isempty(GFPbgMask);
    pause(2)
    figure
    imagesc(mat2gray(GFPimages));
    GFPbgMask=zeros(sx,sy);
    disp('select background for GFP images')
    bf=imfreehand();
    useMask=createMask(bf);
    GFPbgMask=GFPbgMask+useMask;
end
GFPbgimages=GFPimages.*GFPbgMask;
bgIntensity=mean2(GFPbgimages);
DataOut.GFPbgMask = GFPbgMask;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Background Cell for Calcium Image Analysis

if isempty(CabgMask);
    pause(2)
    figure
    imshow(mat2gray(RawImg));
    CabgMask=zeros(sx,sy);
    disp('select background for CALCIUM images')
    bf=imfreehand();
    useMask=createMask(bf);
    CabgMask=CabgMask+useMask;
end

for i=1:sz
    cabgimages(:,:,i)=medfilt2(imagesRaw(:,:,i),[5 5]);
end
DataOut.CabgMask = CabgMask;
thresh = max(max(CabgMask));

%% Calculate the std dev of the background noise%%%%%%%%%%%%
stddevMapbg = zeros(sx,sy);
[x,y]=find(CabgMask==1);
t=1:sz;
for i=1:length(x)
    dat=cabgimages(x(i)-1:x(i)+1,y(i)-1:y(i)+1,:);
    
    dat=mean(dat,2);
    dat=mean(dat,1);
    dat=dat(:);
    
    f=polyfit(t',dat,2);
    feval=polyval(f,t');
    dat2 = dat./feval;
    
    %%try standard deviation
    stddevMapbg(x(i),y(i))=std(dat2);
end

NoiseMap=stddevMapbg.*CabgMask;
StdThresh = max(max(NoiseMap));

%% Calculate Activity %%%%%%%%%%

allactivity=zeros(sx,sy);
for j =1:size(allmasks, 3)
    images = imagesRaw;
    GFPimages = GFPimagesRaw;
    mask=allmasks(:,:,j);
    labeledmask = labeledmask + mask*j;
    %figure; imshow(mask);
    labelbw=bwlabel(mask,4);
    numAreas=max(labelbw(:));
    currArea=1;
    
    %%%%% Filter images for noise %%%%%%%%%
    for i=1:sz
        images(:,:,i)=medfilt2(imagesRaw(:,:,i),[5 5]).*mask;
    end
    GFPimages=GFPimages.*mask;
    
    cellTC(:,j)=shiftdim(mean(mean(images,2),1),1);
    %    f=polyfit(t',dat,2);
    %     feval=polyval(f,t');
    %     dat=dat-feval;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Area=size(nonzeros(mask),1);
    
    %% CALCULATING MAPS
    timeVec=T;
    Fs=mean(diff(timeVec));
    Fs=1/Fs;
    
    Mean = squeeze(mean(mean(images,1),2));
    MeanGFP = mean2(GFPimages);
    
    GFPintensity(j)= mean2(MeanGFP);
    GFPintensityminusbg(j)=GFPintensity(j)-bgIntensity;
    GFPintensityoverbg(j)=GFPintensity(j)/bgIntensity;
    
    oscMap=zeros(sx,sy);
    peakMap=zeros(sx,sy);
    troughMap=zeros(sx,sy);
    timeMat=zeros(sx,sy);
    peakAmpMap=zeros(sx,sy);
    corrMap=nan(sx,sy);
    coeffcorrMap=zeros(sx,sy);
    stddevMap = zeros(sx,sy);
    maxpeakMap = zeros(sx,sy);
    
    [x,y]=find(mask==1);
    
    for i=1:length(x)
        dat=images(x(i)-1:x(i)+1,y(i)-1:y(i)+1,:);
        
        dat=mean(dat,2);
        dat=mean(dat,1);
        dat=dat(:);
        
        f=polyfit(t',dat,2);
        feval=polyval(f,t');
        dat2 = dat./feval;
        dat=dat-feval;
        
        %%try standard deviation
        stddevMap(x(i),y(i))=std(dat2);
        [peakLoc1,peakAmp1]=peakfinder(dat2,(max(dat2)-min(dat2))/5);
        maxpeakMap(x(i),y(i)) = (max(peakAmp1)-1);
        
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
            corrMap(x(i),y(i))=lags(maxCL)/Fs;
            coeffcorrMap(x(i),y(i))= maxCV;
        end
        
    end
    
    DataOut.PeakAmpMap(:,:,j)=peakAmpMap;
    DataOut.maxpeakMap(:,:,j)=maxpeakMap;
    
    %% CALCULATING ACTIVITY and GENERATING ACTIVITY MAPS
    
    %%%%%% thresholding areas above 2*silent cell threshold%%%%%%%%%%%%%
    ActMap=zeros(sx,sy);
    ActMap(logical(mask))=1;
    ActMap(peakAmpMap<2*thresh)=0;
    ActMask=logical(ActMap);
    %% MODIFYING ACTIVITY MAPS
    
    %%%%%% removing regions in the no signal mask, outside of the activity mask,
    %and the image mask%%%%%%%%%%%%
    CoorMap=coeffcorrMap;
    %CoorMap=imfilter(CoorMap, [1 1]);
    CoorMap(~ActMask)=0;
    CoorMap(find(CoorMap<.25))=0;
    
    IntPeakMap = peakAmpMap;
    IntPeakMap(~mask)=0;
    IntPeakMap(IntPeakMap==0)=NaN;
    AveragePeakAmp(j) = nanmean(nanmean(IntPeakMap));
    
    maxpeakMap(~ActMask)=0;
    maxpeakMap = maxpeakMap/StdThresh;
    
    %Map without filtering
    FullActivityMap=zeros(sx,sy);
    FullActivityMap = maxpeakMap;
    FullActivityMap(maxpeakMap==0)=NaN;
    
    allactivity = allactivity+CoorMap;
    
    ActiveArea(j)=size(nonzeros(CoorMap),1);
    NewActiveArea(j)=size(nonzeros(maxpeakMap),1);
    TotalArea(j) = length(nonzeros(mask));
    
    CoorMap(CoorMap==0)=NaN;
    AverageCoor(j) = nanmean(nanmean(CoorMap));
    
    maxpeakMap(find(abs(maxpeakMap)<maxPeakThresh))=0;
    
    %CoorMap(find(CoorMap<activitycutoff))=0;%%%%%%%threshold for activity
    %DataOut.ActivityMap(:,:,i) = CoorMap
end
cellcount = j;%size(cellTC,2)
DataOut.labeledmask = labeledmask;
DataOut.allactivitymap = allactivity;
DataOut.ActiveArea = ActiveArea;
DataOut.AverageCoor = AverageCoor;
DataOut.TotalArea = TotalArea;
DataOut.RatioActive = ActiveArea./(TotalArea);
DataOut.NewRatioActive = NewActiveArea./(TotalArea);
DataOut.GFPintensity = GFPintensity;
DataOut.GFPintensityminusbg = GFPintensityminusbg;
DataOut.GFPintensityoverbg = GFPintensityoverbg;
DataOut.cellTC = cellTC;
DataOut.AveragePeakAmp = AveragePeakAmp;

numofplots = ceil(cellcount/7);
legendstring = cell(1, cellcount);
for i = 1:cellcount
    legendstring{i} = ['cell ' num2str(i)];
end

pause(2)
for kk = 1:numofplots
    if kk == numofplots
        figure
        plot(cellTC(:, 7*kk-6:end))
        legend(legendstring(7*kk-6:end), 'Location','northeastoutside')
    else
        figure
        plot(cellTC(:,7*kk-6:7*kk))
        legend(legendstring(7*kk-6:7*kk), 'Location','northeastoutside')
    end
    title('Timecourses of Selected Cells')
end

%% PLOTTING ACTIVITY MAP
pause(2)
newmap = hsv(cellcount);       %starting map          %how big is it?  %2/3 of way through
newmap = [[0 0 0];newmap];        %set that position to white
colormap(newmap);
figure
%overlay = imoverlay(mat2gray(RawImg), labeledmask, 'colormap', 'hsv(30)', 'facealpha', .75, 'ZeroAlpha', 0);
%imshow(overlay)
%colorbar('Ticks',0:0.5:1,'TickLabels',{'off','on'})
imagesc(labeledmask);
%,'TickLabels',{'off','on'})
%colormap(hsv);
colormap(newmap);
colorbar('Ticks',0:cellcount)
%colorbar('Ticks',0:cellcount)%colorbar
title('Area Map')

allactivity(find(allactivity==0))=NaN;
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
t = linspace(0,1,n)';
r = interp1(t0,colors(:,1),t);
g = interp1(t0,colors(:,2),t);
b = interp1(t0,colors(:,3),t);
cmap = [r,g,b];

pause(2)
imoverlayNurin(mat2gray(RawImg),allactivity,[0, 1],[],cmap,0.3)
colorbar
title('Activity Map')