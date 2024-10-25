function [ dataOut ] = HIanalyzeUNurin( images,timeVec,flipBool,mask )
%HIANALYZE Summary of this function goes here
%   Detailed explanation goes here

%images=cat(3,images{:});
images=double(images);

imagesRaw=mat2gray(images);
%xrealigns images if Islets are moving
%images=registerImages(images);

images=ImgPrep(imagesRaw);

% size(images)

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
%Mask first image
if isempty(mask);
imshow(mat2gray(imagesRaw(:,:,1)));
%numAreas=input('Enter number of areas');
numAreas=1;
mask=zeros(sx,sy);

% level=graythresh(mat2gray(images(:,:,1)));
% mask=im2bw(mat2gray(images(:,:,1)),level);
for i=1:numAreas
bf=imfreehand();
useMask=createMask(bf);
mask=mask+useMask;
end
end
labelbw=bwlabel(mask,4);
numAreas=max(labelbw(:));
currArea=1;

for i=1:sz
   images(:,:,i)=images(:,:,i).*mask;
end


DM=mean(images,1);
DM=mean(DM,2);




for i=1:size(images,3)
    DMean(i)=mean2(images(:,:,i));
end
DF=fft(DMean);
DF(1)=0;
[peakL,peakV]=peakfinder(abs(DF(1:20)));
               
 [mL,mX]=max(peakV);
 
 funFreq=peakL(mX);


freqMap=zeros(sx,sy);
powerMap=zeros(sx,sy);
phaseMap=zeros(sx,sy);
oscMap=zeros(sx,sy);
flipBool=false;
peakMap=zeros(sx,sy);
timeMat=zeros(sx,sy);
peakAmpMap=zeros(sx,sy);
corrMap=nan(sx,sy);
%H=fspecial('average',[3 3]);
%images=imfilter(images,H);
while currArea<=numAreas
    
   [x,y]=find(mask==currArea);
  % [Y,f]=FFT(images,timeVec);
      
       for i=1:length(x)
               
              %
               
               %[maxH]=peakfinder(2*abs(cY(5:length(f))));
               %maxH=maxH+5;
               dat=images(x(i)-1:x(i)+1,y(i)-1:y(i)+1,:);
               
               %dat=images(x(i),y(i),:);
               %dat=reshape(dat,sz,size(dat,1)*size(dat,2));
               dat=mean(dat,2);
               dat=mean(dat,1);
               dat=dat(:);
               %dat=smooth(dat);
               dat=dat(:);
               if flipBool==true
                   dat=1./dat;
               end
               %datFit=detrend(dat);
               %datFitCurve=polyval(datFit,1:length(dat));
               %dat=dat-datFitCurve';
               %dat=smooth(dat);
               dat=detrend(dat);
               Y=fft(dat,NFFT)/L;
               
               
               
               
               F1=2*abs(Y(2:NFFT/2+1));
               [fPeak,fVal]=peakfinder(F1,(max(F1)-min(F1))/4);
               numpeaks=size(fPeak,1);
               [sF,sFX]=sort(F1);
               
               [nm,mm]=max(fVal);
                   
               maxFreq=fPeak(mm);
%                if maxFreq==1
%                    maxFreq=fPeak(2);
%                end
               P=angle(F1(maxFreq));
               maxFreq=TF(maxFreq);
               maxPower=max(fVal);
               freqMap(x(i),y(i))=maxFreq;
%                if (maxFreq==1)
%                    fPeak(1)=[];
%                    fVal(1)=[];
%                    maxFreq=fPeak(find(max(fVal)));
%                    maxPower=max(fVal);
%                    P=angle(F(maxFreq));
%                    
%                end
%                if (maxPower<10*mean(F1(:)))
%                   maxFreq=find(max(F1(:))==F1(:));
%                maxPower=max(F1(:));
%                    P=0;
%                end
     
               
               [peakLoc,peakAmp]=peakfinder(dat,(max(dat)-min(dat))/5);
              % freqMap(x(i),y(i))=length(peakLoc)/max(timeVec);
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
               end
               
               if isempty(peakLoc)
                   numpeak=0;
               end
               oscMap(x(i),y(i))=numpeaks/(max(timeVec));
               
              [c lags]=xcov(dat,DMean,'coeff');
              [maxCV maxCL]=max(c);
              
              phaseMap(x(i),y(i))=angle(Y(funFreq));
              
               if ~isnan(dat(1))
               corrMap(x(i),y(i))=lags(maxCL)/Fs;
               coeffcorrMap(x(i),y(i))= maxCV;
               end
               
          
       end
       
   currArea=currArea+1;
   
end


%% Added by NL to measure % correlation



[M,N]=size(coeffcorrMap);

C25=0;
C50=0;
C75=0;
for m=1:M
    for n=1:N

       
        
        
        if coeffcorrMap(m,n)>=0.25
            C25=C25+1;
        end
        
                if coeffcorrMap(m,n)>=0.5
            C50=C50+1;
                end
        
                        if coeffcorrMap(m,n)>=0.75
            C75=C75+1;
        end
        
        
    end
end

perC25=C25/histcounts(nonzeros(mask),1);
perC50=C50/histcounts(nonzeros(mask),1);
perC75=C75/histcounts(nonzeros(mask),1);


%%

phaseMap=phaseMap/(2*pi*TF(funFreq));

minPhaseMap=min(phaseMap(:));

phaseMap(find(phaseMap==0))=NaN;

phaseMap=phaseMap+abs(minPhaseMap);

dataOut.corrMap=corrMap;
dataOut.coeffcorrMap=coeffcorrMap;
dataOut.freqMap=freqMap;
dataOut.powerMap=powerMap;
dataOut.phaseMap=phaseMap;
dataOut.oscMap=oscMap;
dataOut.timeMat=timeMat;
dataOut.peakMat=peakMap;
dataOut.peakAmpMat=peakAmpMap;
%dataOut.images=images;
dataOut.mask=mask;

dataOut.percentCorr=[perC25 perC50 perC75];



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


[K,L]=size(peakAmpMap);

for k=1:K
    for l=1:L
   if peakAmpMap(k,l)==0
       peakAmpMap(k,l)=NaN;
   end
    end
end
pause(2)

%dataOut=CoordinatedAreas2mM(dataOut,images);


end

