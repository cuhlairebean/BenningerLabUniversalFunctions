function [ dataOut ] = GetPhaseMap(Data, images,timeVec,mask, funFreq )
dataOut = Data;
imageMean=mean(images,3);
img=images(:,:,1);
sx=size(img,1);
sy=size(img,2);
sz=size(images,3);
L=sz;
NFFT=2^nextpow2(L);
Fs=mean(diff(timeVec));
Fs=1/Fs;
TF=Fs./2*linspace(0,1,NFFT/2+1);
TF(1)=[];
TF=1./TF;

DM=mean(images,1);
DM=mean(DM,2);

for i=1:size(images,3)
    DMean(i)=mean2(images(:,:,i));
end
DF=fft(DMean);
DF(1)=0;
[peakL,peakV]=peakfinder(abs(DF(1:18)));

[mL,mX]=max(peakV);
if funFreq == 0
    funFreq=peakL(mX);
end

freqMap=zeros(sx,sy);
powerMap=zeros(sx,sy);
phaseMap=zeros(sx,sy);
oscMap=zeros(sx,sy);
flipBool=false;
peakMap=zeros(sx,sy);
timeMat=zeros(sx,sy);
peakAmpMap=zeros(sx,sy);
corrMap=nan(sx,sy);

[x,y]=find(mask==1);

for i=1:length(x)
    dat=images(x(i)-1:x(i)+1,y(i)-1:y(i)+1,:);
    dat=mean(dat,2);
    dat=mean(dat,1);
    dat=dat(:);
    dat=dat(:);
    
    dat=detrend(dat);
    Y=fft(dat,NFFT)/L;
    
    F1=2*abs(Y(2:NFFT/2+1));
    [fPeak,fVal]=peakfinder(F1,(max(F1)-min(F1))/4);
    numpeaks=size(fPeak,1);
    [sF,sFX]=sort(F1);
    
    [nm,mm]=max(fVal);
    
    maxFreq=fPeak(mm);
    
    P=angle(F1(maxFreq));
    maxFreq=TF(maxFreq);
    maxPower=max(fVal);
    if isempty(maxFreq)
        maxFreq = 0;
    end
    freqMap(x(i),y(i))=maxFreq;
    
    [peakLoc,peakAmp]=peakfinder(dat,(max(dat)-min(dat))/5);
    [trough,troughAmp]=peakfinder(dat,(max(dat)-min(dat))/5,[],-1);
    numpeaks=size(peakLoc,1);
    peakMap(x(i),y(i))=numpeaks;
    
    if length(peakLoc)>1
        timeMat(x(i),y(i))=mean(diff(peakLoc,1));
        peakAmpMap(x(i),y(i))=(mean(peakAmp)-mean(troughAmp));
        phaseMap(x(i),y(i))=mean(peakLoc,1);
    end
    if length(peakLoc)==1
        timeMat(x(i),y(i))=peakLoc;
        peakAmpMap(x(i),y(i))=peakAmp-mean(troughAmp);
        phaseMap(x(i),y(i))=peakLoc;
    end
    if isempty(peakLoc)
        timeMat(x(i),y(i))=0;
        peakAmpMap(x(i),y(i))=0;
        numpeak=0;
    end
    oscMap(x(i),y(i))=numpeaks/(max(timeVec));
    
    [c lags]=xcov(dat,DMean,'coeff');
    [maxCV maxCL]=max(c);
    
    phaseMap(x(i),y(i))=angle(Y(funFreq));
    
    if ~isnan(dat(1))
        corrMap(x(i),y(i))=lags(maxCL)/Fs;
    end
end


phaseMap=phaseMap/(2*pi*TF(funFreq));

minPhaseMap=min(phaseMap(:));

phaseMap(find(phaseMap==0))=NaN;

%phaseMap=phaseMap+abs(minPhaseMap);

dataOut.corrMap=corrMap;
dataOut.freqMap=freqMap;
dataOut.powerMap=powerMap;
dataOut.phaseMap=phaseMap;
dataOut.oscMap=oscMap;
dataOut.timeMat=timeMat;
dataOut.peakMat=peakMap;
dataOut.peakAmpMat=peakAmpMap;
dataOut.mask=mask;
dataOut.funFreq=funFreq;

end

