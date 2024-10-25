function [DataOut]=adjstacksfortime(data,cachannel,howmanychannel,zstacks, starttime, endtime)
%% LOADING IMAGES
addpath('C:\Users\dwuletj\Documents\GitHub\UniversalCode\')
DataOut = data;

%% full correlation map
if starttime == -1
    starttime = 1;
end
if endtime == -1
    endtime = size(DataOut.images,3);
end
T=starttime:endtime;
eachstack =  1.1542;

realTime = 1:eachstack:(endtime*zstacks*eachstack);
for zz=1:zstacks
   timemat(zz,:) = realTime(zz:zstacks:end); 
end

interpolatedimages = zeros(512,512,length(realTime),zstacks);
for zz=1:zstacks
    xq = realTime;
    T = timemat(zz,:);
    for i=1:512
        for j=1:512
            stack = DataOut.images(i,j,starttime:endtime,zz);
            stack = stack(:)';
            vq1 = interp1(T,stack,xq,'linear','extrap');
            interpolatedimages(i,j,:,zz)=vq1;
        end
    end
end
DataOut.interpolatedimages = interpolatedimages;

end