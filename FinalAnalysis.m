function [ dat ] = FinalAnalysis( dat)
%FINALANALYSIS Summary of this function goes here
%   Detailed explanation goes here


BArea=regionprops(dat.Img,'FilledImage');


for i=1:length(BArea)
BoundArea(i)=length(find(BArea(i).FilledImage));
    
end




dat.AMCA=max(dat.Area);
dat.RMCA=dat.AMCA/(length(find(dat.mask)));
dat.AAA=length(find(dat.Img));
dat.RAA=(dat.AAA)/length(find(dat.mask));
dat.period=mean(dat.timeMat(find(dat.Img)));
dat.AvgArea=mean(dat.Area(find(dat.Area)));

AreaUse=dat.Area;
if isempty(AreaUse)
    AreaUse=0; 
end
AreaUse=sort(AreaUse);

fP=length(AreaUse)-ceil(length(AreaUse)*0.05);

if fP~=0
dat.SumFifthArea=sum(AreaUse(fP:end));
dat.SumFifthAreaPerc=dat.SumFifthArea/length(find(dat.mask));
else
dat.SumFifthArea=0;
dat.SumFifthAreaPerc=0; 

end
    
end

