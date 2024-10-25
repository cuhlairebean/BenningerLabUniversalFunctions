function [ Dat ] = FrequencyAnalysis( Dat )
%FREQUENCYANALYSIS Summary of this function goes here
%   Detailed explanation goes here

R=bfopen(Dat.Location);
try
pixelSize=R{4}.getPixelsPhysicalSizeX(0).getValue();
catch
    pixelSize=0.64
end
    pics=R{1};
pics=pics(:,1);



try R{4}.getPlaneDeltaT(0,2)

if ~isempty(R{4}.getPlaneDeltaT(0,0))
for i=1:length(pics)
T(i)=R{4}.getPlaneDeltaT(0,i-1);
end
else
    T=1:size(pics,3);
end

catch

    try
for i=1:length(pics)
T(i)=R{4}.getPlaneDeltaT(i-1,0);
end    
    catch
        
        T=0:0.5:size(pics,3);
    end
end

T=double(T);
Time=T;
 for i=1:length(pics)
IMG(:,:,i)=pics{i};
 end
pics={};
images=double(IMG);
images=mat2gray(images);

IMG=[];

sz=length(Time);

Area=Dat.Area;
peakNum=[];
CAreas=Dat.Img;
%CAreas=repmat(CAreas,1,1,sz);
Fs=mean(diff(Time));
Fs=round(Fs);
Fs=1./Fs;
for i=1:length(Area)
    UPos=find(CAreas==i);
    if ~isempty(UPos)
        for j=1:sz
           Data=images(:,:,j);
           Data=Data(UPos);
           DPlot(i,j)=mean(Data);
        end
      
        %plot(DPlot)
        DPlot(i,:)=detrend(DPlot(i,:));
        %[pxx f]=pwelch(DPlot(i,:),[],[],[],Fs);
        %pxx=pxx./max(pxx);
        %PSectrum(i,:)=pxx;
        DSmooth=DPlot(i,:);
        %for x=1:5
        %    DSmooth=smooth(DSmooth);
        %end
        peaks=peakfinder(DSmooth);
        if ~isempty(peaks)
        if peaks(1)==1
            peaks(1)=[];
        end
        end
        if ~isempty(peaks)
        if peaks(end)==sz;
            peaks(end)=[];
        end
        else
            peaks=[];
        end
        peakNum(i)=length(peaks);  
        
        midData=(max(DPlot(i,:))+min(DPlot(i,:)))/2;
        DCycle(i)=length(find(DPlot(i,:)>midData))/sz;
        
    end
    
    
end

if ~isempty(peakNum)
%PSectrum(find(PSectrum==0))=NaN;
peakNum(find(peakNum==0))=NaN; 
Dat.Frequency=peakNum/Time(end);
DCycle(find(DCycle==0))=NaN;
Dat.DCycle=DCycle;
else
Dat.Frequency=[];
Dat.DCycle=[];
end
%Dat.Spectrum=PSectrum;

end

