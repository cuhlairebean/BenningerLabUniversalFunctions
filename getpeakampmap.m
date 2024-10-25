function [peakAmpMap, stddevMap, corrMap, coeffcorrMap]= getpeakampmap(T, images, mask)%% CALCULATING MAPS
[sx, sy, sz] = size(images);

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
    startx = x(i)-5;
    endx = x(i)+5;
    starty = y(i)-5;
    endy = y(i)+5;
    
    if startx < 1
        startx = 1;
    end
    if endx > 512
        endx = 512;
    end
    if starty < 1
        starty = 1;
    end
    if endy > 512
        endy = 512;
    end
    
    dat = images(startx:endx,starty:endy,:);
    
    dat = mean(dat,2);
    dat = mean(dat,1);
    dat = dat(:);
    
    f = polyfit(t',dat,2);
    feval = polyval(f,t');
    dat2 = dat./feval;
    %dat = dat-feval;
    dat = detrend(dat);
    
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

figure;
imagesc(peakAmpMap);
colorbar;
title('Peak Amp Map');

