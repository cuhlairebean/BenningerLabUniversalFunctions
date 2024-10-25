function [Dat] = OpenFile_JD( data )
% open file with bfopen and remove any unneeded info
if data.newmask == 1
    mask=[];
end
R=bfopen(data.Location);
pics=R{1};
pics=pics(:,1);% remove directory info

%get time slices
try R{4}.getPlaneDeltaT(0,2); 
    if ~isempty(R{4}.getPlaneDeltaT(0,0))
        for i=1:length(pics)
            T(i)=R{4}.getPlaneDeltaT(0,i-1);
        end
    else
        T=1:size(pics,3);
    end
catch
    T=0:0.5:length(pics);    
end

T=double(T);

%change variable for pics
for i=1:length(pics)
    IMG(:,:,i)=pics{i};
end
pics={};
Dat.images = IMG;
Dat.TSize = T;

[Dat]=HIanalyze_JD(IMG,T,0,mask);
Dat.Location=data.Location;
Dat.Image=IMG(:,:,1);
%Dat=CoordinatedAreas(Dat);
end

