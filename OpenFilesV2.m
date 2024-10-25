function [ Dat ] = OpenFiles(data  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if length(fieldnames(data))==1
    data.mask=[];
end
R=bfopen(data.Location);


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
    
% for i=1:length(pics)
% T(i)=R{4}.getPlaneDeltaT(i-1,0);
% end    
    
T=0:0.5:length(pics);    
end

T=double(T);

 for i=1:length(pics)
IMG(:,:,i)=pics{i};
 end
pics={};


mask=data.mask;
[Dat]=HIanalyzeU(IMG,T,0,mask);
Dat.Location=data.Location;
Dat.Image=IMG(:,:,1);
Dat=CoordinatedAreas(Dat);
%Dat=FrequencyAnalysis(Dat);
%aDat=WaveVelocity(Dat);

end

