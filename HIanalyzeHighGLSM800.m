function [ dataOut ] = HIanalyzeHighGLSM800( images,timeVec,flipBool,direc,stimMask,mask )

figure
plot(diff(timeVec))
[stimTime]=peakfinder(diff(timeVec),0,2.1,1);
images=mat2gray(images);
stimTime=stimTime %+1;
%stimTime=[60 70 80 90];
P1=mean(images,3); % average over time
P2=graythresh(P1)/1.7;
%mask=im2bw(P1,P2);
h=fspecial('average');
images=imfilter(images,h);
pulseFreq=mean(diff(stimTime));
%stimTime=[11 21 31];
%images=images(:,:,1:40);
%timeVec=timeVec(1:40);

sx=size(images,1);
sy=size(images,2);
sz=size(images,3);
%mask=stimMask;

if mask==0
    figure
    imshow(images(:,:,1));
    mask=imfreehand();
    mask=createMask(mask);
end
img=mean(images,2);
img=mean(img,1);
img=img(:);
T=mean(diff(timeVec));
Fs=1/T;
L=length(timeVec);
NFFT=2^nextpow2(L);
%pFit=polyfit(timeVec,img,1);
%fitLine=polyval(pFit,timeVec);
%img=img-fitLine;
[imgPeakLocs]=peakfinder(img,(max(img)-min(img))/3);
if isempty(stimTime)
    stimTime=imgPeakLocs;
end
f=Fs/2*linspace(0,1,NFFT/2+1);
YF=fft(img);
[Fun,yFun]=peakfinder(abs(YF(2:NFFT/2+1)));
if length(Fun)>2
    fundFreq=Fun(2);
else
    fundFreq=Fun(1);
end
for i=1:sz
    images(:,:,i)=images(:,:,i).*mask;
end
freqMap=zeros(sx,sy,2);
powerMap=zeros(sx,sy,2);
phaseMap=zeros(sx,sy,2);
oscMap=zeros(sx,sy,2);
corrMap=NaN(sx,sy,1);
lagMap=zeros(sx,sy,2);
ampMap=zeros(sx,sy,2);
pulseMap=nan(sx,sy,2);
drivingMap=zeros(sx,sy);
ratioMap=zeros(sx,sy);
corrMapSaw=NaN(sx,sy,1);
corrMapDrive=NaN(sx,sy,1);
T0=1;
TE=stimTime(1);
powerRatioMap=NaN(sx,sy);
[x,y]=find(mask);
%[Y,f]=FFT(images,timeVec);
for i=1:length(x)
    
    %
    dat=images(x(i)-3:x(i)+3,y(i)-3:y(i)+3,T0:TE);
    
    %dat=images(x(i),y(i),:);
    %dat=reshape(dat,sz,size(dat,1)*size(dat,2));
    dat=mean(dat,2);
    dat=mean(dat,1);
    dat=dat(:);
    %
    %
    YF=fft(dat);
    YF=abs(YF/length(T0:TE));
    YF(1)=0;
    YF(2)=0;
    %
    [maxPower,maxFreq]=max(abs(YF(1:floor(size(YF)/2))));
    [fPeak,fVal]=peakfinder(dat,(max(dat)-min(dat))/4);
    %
    %
    pulseMap(x(i),y(i),1)=length(fPeak);
    %
    powerMap(x(i),y(i),1)=maxPower;
    freqMap(x(i),y(i),1)=(timeVec(TE)-timeVec(T0))/maxFreq;
    
    
    
end

T0=stimTime(1);
TE=sz;
dmean=images(:,:,T0:TE);
dmean=mean(dmean,1);
dmean=mean(dmean,2);
dmean=dmean(:);
DF=fft(dmean);
DF(1)=0;
DF(2)=0;
DT=diff(timeVec);
DT(find(DT==2))=1;
DT=DT(T0-1:TE-1);

imgFT=fft(img(1:T0-1));
imgFT(1)=0;
[~ imgFreq]=max(imgFT(1:30));
dataOut.naturalFreq=imgFreq./(timeVec(T0-1));
dataOut.naturalPeriod=1./dataOut.naturalFreq;
pulseFT=fft(DT);
pulseFT(1)=0;
[~ pulseFreq]=max(pulseFT(1:20));
[xm,ym]=find(stimMask);
for j=1:length(xm)
    stimMean(1:(TE-T0+1),j)=images(xm(j),ym(j),T0:TE);
end
stimMean=mean(stimMean,2);
% stimMean=nanmean(stimMean,1);
% stimMean=nanmean(stimMean,2);
stimMean=stimMean(:);

stimF=fft(stimMean);
stimF(1)=0;
stimF(2)=0;
stimF(3)=0;
[peakLocs,peakVals]=peakfinder(abs(DF(1:floor(size(DF,1)/2))));
[peakValsSort IDX]=sort(peakVals,'descend');

% drivingFreq=peakLocs(find(IDX==2));
% drivingFreq=length(stimTime);

[mm mn]=max(abs(stimF(1:30)));
drivingFreq=mn;

if abs(mn-length(stimTime))>5
    drivingFreq=length(stimTime);
end

if abs(stimF(drivingFreq))<abs(stimF(drivingFreq+1))
    drivingFreq=drivingFreq+1;
end
dataOut.drivingPeriod=(timeVec(TE)-timeVec(T0))/drivingFreq;

%drivingFreq=length(stimTime);
ST=sawtooth(2*pi*length(stimTime)/(TE-T0)*(1:length(DT)));
ST=fliplr(ST);
for i=1:length(x)
    
    
    dat=images(x(i)-3:x(i)+3,y(i)-3:y(i)+3,T0:TE);
    pdat=images(x(i)-3:x(i)+3,y(i)-3:y(i)+3,1:T0-1);
    %dat=images(x(i),y(i),:);
    %dat=reshape(dat,sz,size(dat,1)*size(dat,2));
    dat=mean(dat,2);
    dat=mean(dat,1);
    dat=dat(:);
    %
    pdat=mean(pdat,2);
    pdat=mean(pdat,1);
    pdat=pdat(:);
    
    YP=fft(pdat);
    YP=abs(YP);
    %
    YF=fft(dat);
    YF=abs(YF);
    
    YF(1)=0;
    YF(2)=0;
    %
    [maxPower,maxFreq]=max(abs(YF(1:floor(size(YF,1)/2))));
    [fPeak,fVal]=peakfinder(dat,(max(dat)-min(dat))/5);
    
    
    stimTimeShift=stimTime-(stimTime(1)-1);
    stimTimeShift2=stimTimeShift+1;
    xMember=ismember(stimTimeShift,fPeak);
    xMember2=ismember(stimTimeShift2,fPeak);
    
    
    pulseMap(x(i),y(i),2)=(length(find(xMember))+length(find(xMember2)))/length(stimTimeShift);
    % [mm mn]=max(abs(YF(1:40)));
    %if mn==drivingFreq
    %  pulseMap(x(i),y(i),2)=1;
    %end
    powerMap(x(i),y(i),2)=maxPower;
    [c,lags]=xcov(DT,dat,'coeff');
    corrMap(x(i),y(i),1)=max(c);
    [c,lags]=xcov(ST,dat,'coeff');
    corrMapSaw(x(i),y(i),1)=max(c);
    freqMap(x(i),y(i),2)=(timeVec(TE)-timeVec(T0))/maxFreq;
    drivingMap(x(i),y(i))=abs(YF(drivingFreq));
    [c,lags]=xcov(stimMean,dat,'coeff');
    corrMapDrive(x(i),y(i))=max(c);
    powerRatioMap(x(i),y(i))=YF(pulseFreq)/YP(imgFreq);
    
end

dataOut.drivingMap=drivingMap;
dataOut.powerMap=powerMap;
dataOut.freqMat=freqMap;
dataOut.pulseMap=pulseMap;
dataOut.ratioMap=drivingMap./powerMap(:,:,1);
dataOut.corrMap=corrMap;
dataOut.amplitudeAvg=HighGAnalyze(dataOut);
dataOut.corrArea=length(find(corrMap>0.5))/length(find(~isnan(corrMap)));
dataOut.pulseArea=length(find(pulseMap(:,:,2)>0.7))/length(find(mask));
dataOut.corrMapSaw=corrMapSaw;
dataOut.corrMapDrive=corrMapDrive;
dataOut.powerRatioMap=powerRatioMap;
dataOut.mask=mask;
end


