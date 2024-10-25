function [ data ] = CoordinatedAreas( data )
%COORDINATEDAREAS Summary of this function goes here
%   Detailed explanation goes here

M=data.peakMat;
T=data.timeMat;
Amp=data.peakAmpMat;
Phase=data.phaseMap;
F=data.freqMap;
P=data.peakMat;
Power=data.powerMap;
Power=Power*1000;

MT=M+sqrt(-1)*T;
MT=abs(MT);
MT(find(M<=3))=0;
maxMT(MT==inf)=0;
maxM=max(MT(:));
%T(find(Amp<0.0015))=0;
maxT=max(T(:));
%MT(find(MT<45))=0;
%MT=M;
TP=T+sqrt(-1)*P;
TP=abs(TP);

images=data.images;
images=mat2gray(images);
%images=registerImages(images);
Time = data.T;
sz=length(Time);
%images=images(:,:,size(images)/2:end);
% for i=1:size(images,3)
%     imshow(images(:,:,i))
%     pause(0.01)
% end
%
%
% imshow(Power,colormap(hsv(1.2*max(Power(:)))));
% [x,y]=ginput;
% PCut=Power(floor(y),floor(x));

%PCut=1.8;
%P(find(Amp<0.005))=0; %%%%%%%%%%%%%%%----------------------- NURIN ADV CHANGE----------------
%TP(find(Amp<0.003))=0;

%Threshold image by peak amplitudes - usefull when bad SNR
% GMat=mat2gray(Amp);
% GL=graythresh(GMat)/1.5;
% P(find(GMat<GL))=0;


%T(find(Power<PCut))=0;
V=2:4:max(P(:))+3;
%V=3:1.5:maxM;
MTAll=0*MT;
count=1;
%f=fspecial('average',[2 2]);
f2=fspecial('average',[2 2]);
P=imfilter(P,f2);
if length(V)==1
    V=[V 7];
end

SN=size(T,1)*size(T,2);
sz=size(images,3);
for i=1:length(V)-1
    A=[];
    stats=[];
    MT2=0*MT;
    MT2(find(P>=V(i) & P<V(i+1)))=1;
    MT3=MT2;
    %imshow(MT2)
    MT3=bwlabeln(MT2,8);
    
    MT3=MT3+max(MTAll(:));
    %MT3(find(MT3))=MT3(find(MT3))+max(MTAll(:));
    MT3(find(MT2==0))=0;
    MTAll=MTAll+MT3;
end

MTAll=MTAll.*(MTAll & data.mask);
A0=53421;
for z=1:10
    A=[];
    STATS=regionprops(MTAll,'area');
    
    
    if ~isempty(STATS)
        for j=1:length(STATS)
            stats(j)=STATS(j).Area;
        end
        
        useStats=find(stats>30);
        %useStats=find(stats>5);
        noStats=find(stats<30);
        
        for i=1:length(noStats)
             MTAll(find(MTAll==noStats(i)))=0;
        end
        
        for j=1:length(useStats);
            [x,y]=find(MTAll==useStats(j));
            idx=sub2ind(size(images),x,y);
            for k=1:sz-1
                idx(:,k+1)=idx(:,k)+SN;
            end
            A1=images(idx);
            A1=mean(A1);
            A1=A1(:);
            A(:,j)=A1;
        end
        
        if size(A,2)==1
            break;
        end
        
        if ~isempty(A)
            A=LinMatFit(A);
            [c,lags]=xcov(A,'coeff');
            
            CU=c(size(A,1),:);
            CX=size(A,2);
            pairs=reshape(CU,CX,CX);
            %pairs=triu(pairs);
            
            if size(A,2)==2 && pairs(2,1)>0.7 %% -------------------NURIN ADV CHANGE --------------- the coorelation thresh
                MTAll(find(MTAll==useStats(2)))=useStats(1);
                MT3=MTAll;
                break;
            end
            
            
            CT=zeros(CX,CX);
            for k=1:size(pairs,1)
                UP=find(pairs(k,:)>0.7 & pairs(k,:)<0.9999999999);%%--------------NURIN ADV CHANGE-------- 0.85 from 0.7 to 0.85
                if ~isempty(UP)
                    for j=1:length(UP)
                        if length(UP)==1
                            CT(k,UP(j))=UP(j);
                        else
                            
                            CT(k,j)=UP(j);
                        end
                    end
                end
            end
            for k=1:size(CT,1)
                numCoord=length(find(CT(k,:)));
            end
            
            if size(A,2)==A0
                for k=1:size(CT,1)
                    useNums=CT(k,:);
                    useNums=useNums(find(useNums));
                    if ~isempty(useNums)
                        useNums=[k useNums];
                        for j=1:length(useNums)
                            MTAll(find(MTAll==useStats(useNums(j))))=MTMax;
                            CT(find(CT==useNums(j) | CT==j))=0;
                        end
                        MTMax=MTMax+1;
                    end
                end
                
                break;
            else
                A0=size(A,2);
                
            end
            
            MTMax=max(MTAll(:))+1;
            for k=1:size(CT,1)
                useNums=CT(k,:);
                useNums=useNums(find(useNums));
                if ~isempty(useNums)
                    for j=1:length(useNums)
                        MTAll(find(MTAll==useStats(useNums(j))))=MTMax;
                        CT(find(CT==useNums(j) | CT==j))=0;
                    end
                    MTMax=MTMax+1;

                end 
            end
            MT3=MTAll;
        end
    end 
end


A=[];
STATS=regionprops(MTAll,'area');
MT3=MTAll;
MTAll=0*MTAll;
count=1;
for i=1:length(STATS)
    AVal=STATS(i).Area;
    if AVal>20
        MTAll(find(MT3==i))=count;
        count=count+1;
        
    end
    
end

for i=1:max(MTAll(:))
    if mean(Amp(find(MTAll==i)))<0.0007
        MTAll(find(MTAll==i))=0;
    end
end
MTAll2=MTAll;
for i=1:max(MTAll(:))
    MTU=MTAll;
    MTU(find(MTU~=i))=0;
    MTU=imfill(MTU);
    MTU(find(MTU))=i;
    MTAll(find(MTU))=i;
    MTAll2(find(MTAll==i & MTAll2==0))=i;
    %MTAll2(find(MTAll==i))=i;
end

MTAll=MTAll2;
MTAll2(find(P<0.01))=0;
data.Img=MTAll;  %coordinated areas map
A1=regionprops(MTAll,'Area');
A=[];




for i=1:length(A1)
    A(i)=A1(i).Area;
end


data.Area=A;
%data.pixelSize=pixelSize;
%data.metaData=R{4};
%figure
data=FinalAnalysis(data);
for i=1:max(MTAll(:))
    A1=find(MTAll==i);
    for j=1:size(images,3);
        useA=images(:,:,j);
        useADat=useA(A1);
        useADat=mean(useADat(:));
        tCourse(j,i)=useADat;
    end
end

for i=1:max(MTAll(:))
    tCourseH(:,i)=unwrap(angle(hilbert(detrend(tCourse(:,i)))));
end
%N=length(find(~isnan(tCourse(i,:))));
%Z=0;
%Z=1/N*nansum(exp(sqrt(-1)*tCourseH),2);
%data.Z=nan;
%data.AvgZ=mean(abs(Z));
%imshow(MTAll,colormap(hsv(length(A))));
data.Image=images(:,:,1);
%figure
%imshow(MTAll,colormap(hsv(length(A)+10)));
data.AvgFreq=mean(data.freqMap(find(data.freqMap)));

data=FrequencyAnalysis(data);


end

