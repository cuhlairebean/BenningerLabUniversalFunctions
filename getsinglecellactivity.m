function data = getsinglecellactivity(data)
%Get activity for single cell analysis;
images = data.images;
mask = data.mask;
bgmask = data.bgmask;
numCells = max(max(mask));
data.numCells = numCells;
[sx, sy, sz] = size(images);

singlecellmask = zeros(size(mask));
bgcellmask = zeros(size(mask));
ActMask=zeros(size(mask));

timeVec=data.T;
Fs=mean(diff(timeVec));
Fs=1/Fs;

Mean = squeeze(mean(mean(images,1),2));

peakMap=zeros(1,numCells+1);
troughMap=zeros(1,numCells+1);
timeMat=zeros(1,numCells+1);
peakAmpMap=zeros(1,numCells+1);

% corrMap=nan(sx,sy);
% coeffcorrMap=zeros(sx,sy);

t=1:sz;

for i =1:numCells+1
    singlecellmask = zeros(size(mask));
    if i <=numCells
        singlecellmask(find(mask == i))=1;
        figure
        imagesc(singlecellmask)
        title('which cell')
    else
        singlecellmask(find(bgmask == 1))=1;
        figure
        imagesc(singlecellmask)
        title('background cell')
    end
    newimages = images.*singlecellmask;
    newimages(~logical(newimages))=NaN;
    dat = nanmean(newimages,2);
    dat = nanmean(dat,1);
    dat = dat(:);
    data.dat(i,:) = dat;
    figure
    subplot(2,1,1)
    plot(dat)
    if i <=numCells
        title(['Cell ' num2str(i) ' orig data'])
    else
        title(['Background Cell orig data'])
    end
    
    %%untrend and normalize data
%     f = polyfit(t',dat,2);
%     feval = polyval(f,t');
%     data.polyval(i,:) = feval;
%     dat2 = dat./feval;
%     dat2 = (detrend(dat2));
    %     dat3 = (dat-feval);
    dat = (detrend(dat))/mean(dat);
    % display untrended data
    data.dat_adj(i,:) = dat;
    subplot(2,1,2)
    plot(dat)
    title(['adj data'])
    
%     figure
%     plot(dat2)
%     title(['test data'])
    
    %%smooth data and get average frequency
    %output = smoothdata(dat,'movmean',5);
    [freqvecmean, ~] = findavgfreq(dat);
    data.freq(i)=freqvecmean;
    
    %[~,peakAmp1]=peakfinder(dat,(max(dat)-min(dat))/5);
    
    [peakLoc,peakAmp]=peakfinder(dat,(max(dat)-min(dat))/5);
    [trough,troughAmp]=peakfinder(dat,(max(dat)-min(dat))/5,[],-1);
    numpeaks=size(peakLoc,1);
    peakMap(i)=numpeaks;
    troughMap(i)=size(trough,1);
    if length(peakLoc)>1
        timeMat(i)=mean(diff(peakLoc,1));
        peakAmpMap(i)=(mean(peakAmp)-mean(troughAmp));
    end
    if length(peakLoc)==1
        timeMat(i)=peakLoc;
        peakAmpMap(i)=peakAmp-mean(troughAmp);
    end
    if isempty(peakLoc)
        timeMat(i)=0;
        peakAmpMap(i)=0;
    end
end
data.peakAmpMap = peakAmpMap;

thresh = peakAmpMap(end);
data.thresh = thresh;

%% which cells are active?
ActMap=zeros(1,numCells);
%Determine the active area using the bachground threshhold
ActMap(peakAmpMap>1.1*thresh)=1;
data.activity = ActMap;

activecells = find(ActMap);
data.activecells = activecells;
data.notactivecells = find(~ActMap);
data.numActive = length(activecells);
data.numNotActive = length(data.notactivecells);

for i = 1:length(activecells)
    ActMask(find(mask==activecells(i))) = 1;
end
figure
imagesc(ActMask);
title('Active Cells')
data.ActMask = ActMask;

%%getDutyCycle
for i =1:numCells
    tc = data.dat_adj(i,:);
    abovethresh = find(tc>1.1*thresh);
    dc = length(abovethresh)/length(tc);
    data.dutycycle(i) = dc;
    figure
    plot(tc)
    hold on
    plot(thresh*ones(1,length(tc)))
    title(['Cell ' num2str(i) 'duty cycle'])
end

end

