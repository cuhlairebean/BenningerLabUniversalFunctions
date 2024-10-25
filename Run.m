function [DataOut]=Run(data)
%% LOADING IMAGES

if length(fieldnames(data))==1
    mask=[];
    silentcell=[];
else
    mask=data.mask;
    silentcell=data.maskSC;
    thresh=data.thresh;
end

R=bfopen(data.Location);

if length(R(:,4))>1
    pics=R(:,1);
    for i=1:length(pics)
        pics1 = pics{i};
        pics2{i} = pics1{1};
    end
    pics=pics2';
else
    pics=R{1};
    pics=pics(:,1);
end

try R{4}.getPlaneDeltaT(0,2)
    
    if ~isempty(R{4}.getPlaneDeltaT(0,0))
        for i=1:length(pics)
            T(i)=R{4}.getPl/aneDeltaT(0,i-1);
        end
    else
        T=1:size(pics,3);
    end
    
catch
    
    % for i=1:length(pics)
    % T(i)=R{4}.getPlaneDeltaT(i-1,0);
    % end
    
    T=0:1:length(pics);
end

T=double(T);

for i=1:length(pics)
    IMG(:,:,i)=pics{i};
end
pics={};

Images=double(IMG);
RawImg=Images(:,:,1);


%% CODE TO MODIFY TCs
st=400;
%ed=size(Images,3);
ed=800;
T=st:1:ed;
%
images=Images(:,:,st:ed);
DataOut.images=images;
% % % %

%% DECLARING IMAGE PROPERTIES

sx=size(images,1);
sy=size(images,2);
sz=size(images,3);

%% MASKING IMAGES

if isempty(mask);
    figure
    imshow(mat2gray(RawImg));
    %numAreas=input('Enter number of areas');
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
end
DataOut.mask = mask;
labelbw=bwlabel(mask,4);
numAreas=max(labelbw(:));
currArea=1;
imagesRaw=images;
%%%%%
h = fspecial('average', 5);
for i=1:sz
    %    images(:,:,i)=imfilter(images(:,:,i),h,'replicate').*mask;
    %images(:,:,i)=images(:,:,i).*mask;
    images(:,:,i)=medfilt2(images(:,:,i),[5 5]);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% CALCULATING MAPS
timeVec=T;
Fs=mean(diff(timeVec));
Fs=1/Fs;

for i=1:sz
    Mean(i)=mean2(images(:,:,i));
end

peakMap=zeros(sx,sy);
troughMap=zeros(sx,sy);
timeMat=zeros(sx,sy);
peakAmpMap=zeros(sx,sy);
corrMap=nan(sx,sy);
coeffcorrMap=zeros(sx,sy);
t=1:sz;

while currArea<=numAreas
    
    [x,y]=find(mask==currArea);
    
    for i=1:length(x)
        
        
        dat=images(x(i)-1:x(i)+1,y(i)-1:y(i)+1,:);
        
        dat=mean(dat,2);
        dat=mean(dat,1);
        dat=dat(:);
        %                size(dat)
        %                size(t)
        f=polyfit(t',dat,2);
        feval=polyval(f,t');
        
        
        dat=dat-feval;
        
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
        
        
        [c lags]=xcov(dat,Mean','coeff');
        [maxCV maxCL]=max(c);
        
        if ~isnan(dat(1))
            corrMap(x(i),y(i))=lags(maxCL)/Fs;
            coeffcorrMap(x(i),y(i))= maxCV;
        end
        
    end
    currArea=currArea+1;
end
DataOut.xInfo=length(x);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%ACTIVITY CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% APPLYING "SILENT CELL" THRESHOLD
DataOut.PeakAmpMap=peakAmpMap;

%%%%%%Selecting "silent cell" from peak amp map%%%%%%%%%%
if selectBackgroundBoolean
    if isempty(silentcell);
        figure;
        imagesc(peakAmpMap)
        disp('select background cell')
        silentcell=createMask(imfreehand);
        
        SCMap=peakAmpMap.*silentcell;
        mxScMap=max(max(SCMap));
        thresh = mxScMap;
    end
    
    DataOut.maskSC=silentcell;
    DataOut.thresh=thresh;
    %%%% using the silent cell threshold to remove areas with NO signal %%%%%%
    [NoSignalMask, NoSignalArea]=RemovingAreaswNoSignal(peakAmpMap,thresh);
   
else
    if isempty(silentcell);
        figure;
        imshow(mat2gray(RawImg));
        disp('select background for activity threshhold')
        silentcell=createMask(imfreehand);
    end
    
    [x1,y1]=find(silentcell==1);
    backgroundMap = zeros(sx, sy);
    for i=1:length(x1)
        dat=imagesRaw(x1(i)-1:x1(i)+1,y1(i)-1:y1(i)+1,:);
        dat=mean(dat,2);
        dat=mean(dat,1);
        dat=dat(:);
        f=polyfit(t',dat,2);
        feval=polyval(f,t');
        dat=dat-feval;
        [peakLoc,peakAmp]=peakfinder(dat,(max(dat)-min(dat))/2);
        [trough,troughAmp]=peakfinder(dat,(max(dat)-min(dat))/2,[],-1);
        numpeaks=size(peakLoc,1);
        if length(peakLoc)>1
            backgroundMap(x1(i),y1(i))=(mean(peakAmp)-mean(troughAmp));
        end
        if length(peakLoc)==1
            backgroundMap(x1(i),y1(i))=peakAmp-mean(troughAmp);
        end
        if isempty(peakLoc)
            backgroundMap(x1(i),y1(i))=0;
        end
    end
    
    BGMap=backgroundMap.*silentcell;
    mxScMap=max(max(BGMap));
    thresh = mxScMap;
    NoSignalMask = zeros(sx, sy);
    NoSignalMask=logical(NoSignalMask);
    NoSignalArea = 0;
    DataOut.thresh=thresh;
    DataOut.maskSC=[];
    DataOut.silentcell = silentcell;
end

Area=size(nonzeros(mask),1)-NoSignalArea;

% hold on
% dataFit=polyfit(1:sz,threshTC,2);
% dataFit2=polyval(dataFit,1:sz);
% normthreshTC=threshTC-dataFit2;
% plot(normthreshTC)
% legend ('rawTC','normTC')

% peakAmpMap=peakAmpMap./mean2(nonzeros(peakAmpMap));

%% CALCULATING ACTIVITY and GENERATING ACTIVITY MAPS

%%%%%% thresholding areas above 2*silent cell threshold%%%%%%%%%%%%%
ActMap=zeros(sx,sy);
ActMap(find(mask))=1;
ActMap(NoSignalMask)=0;
ActMap(find(peakAmpMap<1.5*thresh))=0;%%%%%%%threshhold for activity
ActMask=logical(ActMap);
ActMap=bwlabeln(im2bw(ActMap));
STATS=regionprops(ActMap,'area');

% DataOut.STATS=STATS;
%    for j=1:length(STATS)
%       stats(j)=STATS(j).Area;
%    end
%
% %    DataOut.stats=stats;
%
% nostats=find(stats<50);
% yesstats=find(stats>50);
% for i=1:length(nostats)
%     ActMap(find(ActMap==nostats(i)))=0;
% end


%% MODIFYING ACTIVITY MAPS

colors = [
    1 0 1 % First element = purple
   0 0 1 % blue
 %  0 1 1 % cyan
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

% CoorArea=ActMask*CorrMap;

CoorMap=zeros(sx,sy);
lCoorMap=zeros(sx,sy);

%%%%%% removing regions in the no signal mask, outside of the activity mask,
%and the image mask%%%%%%%%%%%%
CoorMap=coeffcorrMap;
CoorMap(~ActMask)=0;
CoorMap(find(CoorMap<0.1))=0;%%%%%%%addtional noise

%% CALCULATING DUTY CYCLE

t=st:ed;
tint=st:0.01:ed;

% CoorMapBW=im2bw(CoorMap);
% CoorMapBWl=bwlabel(CoorMapBW,8);
% clear STATS stats nostats yesstats val loc
% STATS=regionprops(CoorMapBWl,'area');
% if length(STATS)>0
%     
%     for k=1:length(STATS)
%         stats(k)=STATS(k).Area;
%     end
%     nostats=find(stats<100);
%     yesstats=find(stats>100);
%     
%     %    CoorMap_rev=CoorMapBWl;
%     
%     if length(yesstats)>0
%         
%         
%         %figure;
%         %title('Duty Cycle')
%         %hold on;
%         for i=1:length(yesstats)
%             
%             
%             maskDC=zeros(sx,sy);
%             
%             maskDC(find(CoorMapBWl==yesstats(i)))=1;
%             
%             mapsDC(:,:,i)=maskDC(:,:);
%             
%             maskDC=repmat(maskDC,[1 1 sz]);
%             Ims=maskDC.*images;
%             Ca=shiftdim(mean(mean(Ims(:,:,:),2),1),1);
%             
%             t=st:ed;
%             Ca=medfilt1(Ca);
%             %             Ca=interp1(t,Ca,tint);
%             fitCa=polyfit(t,Ca,3);
%             corrCa=polyval(fitCa, t);
%             Ca=Ca./corrCa;
%             %           plot(Ca)
%             
%             % drawnow
%             [pLoc, pAmp]=peakfinder(Ca);
%             [tLoc,tAmp]=peakfinder(Ca,[],[],-1);
%             amp=mean(pAmp)-mean(tAmp);
%             %            scatter(pLoc,pAmp)
%             %            scatter(tLoc,tAmp)
%             %            drawnow
%             
%             %
%             
%             time=0;
%             count=1;
%             clear loc val
%             
%             for k=1:length(t)
%                 %     amp=1.0;
%                 amp=(mean(pAmp)+mean(tAmp))/2;
%                 %     if Ca(k)>=amp*0.05+mean(tAmp)
%                 if Ca(k)>=amp
%                     time=time+1;
%                     loc(count)=k;
%                     val(count)=Ca(k);
%                     count=count+1;
%                 end
%             end
%             %scatter(loc,val, 'filled')
%             %drawnow
%             
%             
%             PerOn(i)=time/sz;
%             
%             % scatter(loc,val,'filled')
%             % drawnow
%             % pause(2)
%             
%             
%         end
%         
%         dutycycle=mean(PerOn);
%         
%         % dutymap=DutyCalculator(images, CoorMapOver20_f)
%         %DataOut.DutyMaps=mapsDC;
%     else
%         PerOn=0;
%         dutycycle=0;
%     end
% else
%     PerOn=0;
%     dutycycle=0;
%     
% end
%% PLOTTING ACTIVITY MAP
%h = fspecial('average', 2);
%CoorMap=imfilter(round(CoorMap,1),h,'replicate');
ActiveArea=size(nonzeros(CoorMap),1);
% NoNANCorrMap=CoorMap;
CoorMap(find(CoorMap==0))=NaN;
pause(2)
%imoverlayNurin(mat2gray(imagesRaw(:,:,1)),medfilt2(CoorMap,[10, 10], 'symmetric'),[0, 1],[],cmap,0.3)
imoverlayNurin(mat2gray(imagesRaw(:,:,1)),CoorMap,[0, 1],[],cmap,0.3)
colorbar
%caxis([0 1]) % forces activity map to be scaled between [0-1] - not
%physical because activity has arbitrary units and can be >1
title('Activity Map')
CorrMap(~ActMask)=NaN;


ActMap(~mask)=NaN;
ActMap(NoSignalMask)=NaN;
CorrMap(NoSignalMask)=NaN;

%
% pause(2)
% imoverlay(mat2gray(imagesRaw(:,:,1)),ActMap,[],[],cmap,0.3)
% colorbar

%  h = fspecial('average', 2);
% pause(2)
% imoverlay(mat2gray(imagesRaw(:,:,1)),imfilter(round(CorrMap,1),h,'replicate'),[0, 1],[],cmap,0.3);
% colorbar

% DataOut.round=round(CorrMap,1);

% DataOut.RawImages=imagesRaw;


%% CALCULATING PERCENT ACTIVE AREA

DataOut.RatioActive=ActiveArea/(Area);

%% SAVING ACTIVITY VARIABLES

DataOut.ActMap=ActMap;
DataOut.ActMapValues = CoorMap;
%DataOut.CorrMap=coeffcorrMap;
DataOut.mask=mask;
DataOut.Location=data.Location;
DataOut.NoSignalMask=NoSignalMask;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% *%%%%%%%%%%%%%%%%%%%%%%%%CORRELATION CALCULATIONS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% SELECTING THREE CELLS AND SAVING TIME COURSES

sz_2=ed-st+1;

peakAmpMap(~ActMask)=0;
if length(fieldnames(data))==1
    
    figure;
    imagesc(peakAmpMap)
    disp('select cell 1 for correlation analysis')
    cell1=createMask(imfreehand);
    maskstack = repmat(cell1,[1 1 sz_2]);
    maskedStack = images .* maskstack;
    cell1TC=shiftdim(mean(mean(maskedStack,2),1),1);
    f=polyfit(t,cell1TC,2);
    cell1TC=cell1TC./polyval(f,t);
    %cell1TC=medfilt1(cell1TC,10);
    %
    disp('select cell 2 for correlation analysis')
    cell2=createMask(imfreehand);
    maskstack = repmat(cell2,[1 1 sz_2]);
    maskedStack = images .* maskstack;
    cell2TC=shiftdim(mean(mean(maskedStack,2),1),1);
    f=polyfit(t,cell2TC,2);
    cell2TC=cell2TC./polyval(f,t);
    %cell2TC=medfilt1(cell2TC,10);
    
    disp('select cell 3 for correlation analysis')
    cell3=createMask(imfreehand);
    maskstack = repmat(cell3,[1 1 sz_2]);
    maskedStack = images .* maskstack;
    cell3TC=shiftdim(mean(mean(maskedStack,2),1),1);
    f=polyfit(t,cell3TC,2);
    cell3TC=cell3TC./polyval(f,t);
    %cell3TC=medfilt1(cell3TC,10);
    
  
    figure
    plot(cell1TC,'-b')
    hold on
    plot(cell2TC,'-r')
    plot(cell3TC,'-g')
    legend('cell 1','cell 2','cell 3')
    title('timecourses of selected cells')
    hold off
    
else
    
    CorrCellTCs=data.CorrCellTCs;
    DataOut.CorrCellTCs_i=CorrCellTCs;
    cell1TC=CorrCellTCs(1,:);
    cell2TC=CorrCellTCs(2,:);
    cell3TC=CorrCellTCs(3,:);
    
end

%% GENERATING CORRELATION MAPS

% labelbw=bwlabel(mask,4);
% numAreas=max(labelbw(:));
% currArea=1;

clear STATS stats nostats yesstats val loc


[sx,sy,sz]=size(images);

CorrMap1=zeros(sx,sy);
CorrMap2=zeros(sx,sy);
CorrMap3=zeros(sx,sy);

t=1:sz;

% intermediateMap = medfilt2(peakAmpMap,[10,10]);
% intermediateMap(~mask) = 0;
% intermediateMap(NoSignalMask) = 0;
%
% CoorMapBW=(intermediateMap./mean(nonzeros(intermediateMap)));
% CoorMapBW(~mask) = NaN;
% CoorMapBW(NoSignalMask) = NaN;

LogicalMap=zeros(sx,sy);
% LogicalMap(find(mask==1))=1;
% LogicalMap(NoSignalMask)=0;
LogicalMap(ActMask) = 1; %only run correlation on active area

maskCorr=zeros(sx,sy);
maskCorr(find(LogicalMap==1))=1;
maskCorr=repmat(maskCorr,[1 1 sz]);
Ims=maskCorr.*images;
%
% dat=shiftdim(mean(mean(Ims(:,:,:),2),1),1);
% f=polyfit(t,dat,2);
% dat=dat./polyval(f,t);
% dat=medfilt1(dat,10);
%

[x,y]=find(LogicalMap==1);
t=1:sz;
for i=1:length(x)
    x1=x(i)-5;
    x2=x(i)+5;
    y1=y(i)-5;
    y2=y(i)+5;
    if x1 < 1
        x1=1;
    end
    if x2 > sx
        x2=sx;
    end
    
    if y1 < 1
        y1=1;
    end
    if y2 > sy
        y2=sy;
    end
    dat=Ims(x1:x2,y1:y2,:);
    
    dat=nanmean(dat,2);
    dat=nanmean(dat,1);
    dat=dat(:);
    f=polyfit(t',dat,2);
    feval=polyval(f,t');
    %dat=dat-feval;
    %dat=detrend(dat)';
    dat=dat./polyval(f,t');
    %dat=medfilt1(dat,10);
    dat=dat';
    [c lags]=xcov(dat,cell1TC, 'coeff');
    [maxCV maxCL]=max(c);
    if ~isnan(dat(1))
        % corrMap(x(i),y(i))=lags(maxCL)/Fs;
        CorrMap1(x(i),y(i))= maxCV;
        CorrMap1_2(x(i),y(i))= lags(maxCL);
    end
    
    [c lags]=xcov(dat,cell2TC, 'coeff');
    [maxCV maxCL]=max(c);
    if ~isnan(dat(1))
        CorrMap2(x(i),y(i))= maxCV;
        CorrMap2_2(x(i),y(i))= lags(maxCL);
    end
    
    [c lags]=xcov(dat,cell3TC, 'coeff');
    [maxCV maxCL]=max(c);
    if ~isnan(dat(1))
        CorrMap3(x(i),y(i))= maxCV;
        CorrMap3_2(x(i),y(i))= lags(maxCL);
    end
    
    
end

%%

%CorrMap1(NoSignalMask)=0; CorrMap2(NoSignalMask)=0; CorrMap3(NoSignalMask)=0;

% CorrMap1=imfilter(round(CorrMap1,1),h,'replicate');
% CorrMap2=imfilter(round(CorrMap2,1),h,'replicate');
% CorrMap3=imfilter(round(CorrMap3,1),h,'replicate');
CorrMap1(~LogicalMap)=NaN; CorrMap2(~LogicalMap)=NaN; CorrMap3(~LogicalMap)=NaN;

pause(2)
% figure('rend','painters','pos',[500 10 400 1000])
% subplot(3, 1, 1); imagesc(CorrMap1);colorbar('Ticks',[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1]); title('Correlation Map for Cell 1');
% subplot(3, 1, 2); imagesc(CorrMap2);colorbar; title('Correlation Map for Cell 1');
% subplot(3, 1, 3); imagesc(CorrMap3);colorbar; title('Correlation Map for Cell 1');

DataOut.CorrCellTCs=[cell1TC; cell2TC; cell3TC];
DataOut.CorrCellMaps(:,:,1)=CorrMap1;
DataOut.CorrCellMaps(:,:,2)=CorrMap2;
DataOut.CorrCellMaps(:,:,3)=CorrMap3;

DataOut.CorrCellMaps_2(:,:,1)=CorrMap1_2;
DataOut.CorrCellMaps_2(:,:,2)=CorrMap2_2;
DataOut.CorrCellMaps_2(:,:,3)=CorrMap3_2;



%%


% DataOut.Vals=Vals;

% FIGURES COMMENTED OUT

% pause(2)
% imoverlay(mat2gray(RawImg),CorrMap1,[0,max(max(CorrMap1))],[],cmap,0.3);
% colorbar
%
% pause(2)
% imoverlay(mat2gray(RawImg),CorrMap2,[0,max(max(CorrMap2))],[],cmap,0.3);
% colorbar
%
% pause(2)
% imoverlay(mat2gray(RawImg),CorrMap3,[0,max(max(CorrMap3))],[],cmap,0.3);
% colorbar

%% THRESHOLDING DOWN LOW-CORELLATED REGIONS
%
% CorrMap1(find(CorrMap1<0.1))=0;
ACorr1=size(nonzeros(CorrMap1.*~isnan(CorrMap1)),1);

% CorrMap2(find(CorrMap2<0.1))=0;
ACorr2=size(nonzeros(CorrMap2.*~isnan(CorrMap2)),1);
%
% CorrMap3(find(CorrMap3<0.1))=0;
ACorr3=size(nonzeros(CorrMap3.*~isnan(CorrMap3)),1);

avgCorr1 = nanmean(nanmean(CorrMap1));
%avgCorr1_nozeros = nanmean(nonzeros(CorrMap1));
%sizeCorr1 = size(~isnan(CorrMap1));
avgCorr2 = nanmean(nanmean(CorrMap2));
%avgCorr2_nozeros = nanmean(nonzeros(CorrMap2));
%sizeCorr2 = size(~isnan(CorrMap2));
avgCorr3 = nanmean(nanmean(CorrMap3));
%avgCorr3_nozeros = nanmean(nonzeros(CorrMap3));
%sizeCorr3 = size(~isnan(CorrMap3));

DataOut.AvgCorr = [avgCorr1 avgCorr2 avgCorr3];
%DataOut.AvgCorr_nozeros = [avgCorr1_nozeros avgCorr2_nozeros avgCorr3_nozeros];

%% CALCULATING CORRELATIONS
%DataOut.RawAreaCorr=[ACorr1; ACorr2; ACorr3]./(Area);
%[ACorr,ind]=sort([ACorr1; ACorr2; ACorr3],'descend');
%ACorr=ACorr(1);


names=['CorrMap1'; 'CorrMap2'; 'CorrMap3'];
%CorrMatrixName=char(names(ind(1),:));
%DataOut.HiCorrMatrixNam=CorrMatrixName;

%HiCorrMatrix=eval(CorrMatrixName);
% pause(2)
% imoverlayNurin(mat2gray(RawImg),HiCorrMatrix,[0,1],[],cmap,0.3)
% colorbar
% title 'Correlation Map'

pause(2)
imoverlayNurin(mat2gray(RawImg),CorrMap1,[0,1],[],cmap,0.3)
colorbar
title 'Correlation Map for Cell 1'
pause(2)
imoverlayNurin(mat2gray(RawImg),CorrMap2,[0,1],[],cmap,0.3)
colorbar
title 'Correlation Map for Cell 2'
pause(2)
imoverlayNurin(mat2gray(RawImg),CorrMap3,[0,1],[],cmap,0.3)
colorbar
title 'Correlation Map for Cell 3'
pause(2)
To30_1=zeros(sx,sy);
To50_1=zeros(sx,sy);
To75_1=zeros(sx,sy);

To30_2=zeros(sx,sy);
To50_2=zeros(sx,sy);
To75_2=zeros(sx,sy);

To30_3=zeros(sx,sy);
To50_3=zeros(sx,sy);
To75_3=zeros(sx,sy);

To30_1(find(CorrMap1<0.3 & CorrMap1>0))=0.3;
sTo30(1)=size(nonzeros(To30_1),1);
To50_1(find(CorrMap1<0.7 & CorrMap1>=0.3))=0.5;
sTo50(1)=size(nonzeros(To50_1),1);
To70_1(find(CorrMap1<=1.0 & CorrMap1>=0.7))=1;
sTo70(1)=size(nonzeros(To70_1),1);

To30_2(find(CorrMap2<0.3 & CorrMap2>0))=0.3;
sTo30(2)=size(nonzeros(To30_2),1);
To50_2(find(CorrMap2<0.7 & CorrMap2>=0.3))=0.5;
sTo50(2)=size(nonzeros(To50_2),1);
To70_2(find(CorrMap2<=1.0 & CorrMap2>=0.7))=1;
sTo70(2)=size(nonzeros(To70_2),1);

To30_3(find(CorrMap3<0.3 & CorrMap3>0))=0.3;
sTo30(3)=size(nonzeros(To30_3),1);
To50_3(find(CorrMap3<0.7 & CorrMap3>=0.3))=0.5;
sTo50(3)=size(nonzeros(To50_3),1);
To70_3(find(CorrMap3<=1.0 & CorrMap3>=0.7))=1;
sTo70(3)=size(nonzeros(To70_3),1);


%ACorr=ACorr(1);
%names=['To30'; 'To50'; 'To70'];
%CorrMatrixName=char(names(ind(1),:));
%DataOut.HiCorrRangeMatrixName=CorrMatrixName;
ACorrs=[sTo30; sTo50; sTo70];
ACorrs=ACorrs/(ActiveArea);
DataOut.ACorrs=ACorrs;

%[ACorr,ind]=sort([sTo30; sTo50; sTo70],'descend');
DataOut.RatioCorr_1=(sTo70(1))./(ActiveArea);
DataOut.RatioCorr_2=(sTo70(2))./(ActiveArea);
DataOut.RatioCorr_3=(sTo70(3))./(ActiveArea);

% ind=ind(1);

% RatioCorr=ACorr./(Area);
% DataOut.RatioCorr=RatioCorr;



% pause(2)
% imoverlay(mat2gray(RawImg),eval(CorrMatrixName),[0,1],[],cmap,0.3);
% colorbar


disp('End of analysis.')

end