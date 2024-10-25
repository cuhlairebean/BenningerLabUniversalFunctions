function SmallAreaMap = RemovingSmallAreas( map, areasize )
intmap=bwlabeln(imbinarize(map));

figure
imagesc(imbinarize(intmap))
title('Map of Areas')

STATS=regionprops(intmap,'area');
DataOut.STATS=STATS;
for j=1:length(STATS)
    stats(j)=STATS(j).Area;
end
SmallAreaMap=zeros(size(map));
nostats=find(stats<areasize);
for i=1:length(nostats)
    SmallAreaMap(intmap==nostats(i))=1;
end
    SmallAreaMap = logical(SmallAreaMap);
pause(2)
% figure
% imagesc(SmallAreaMap);
% colormap('jet')
% title('Small Area Map');
end

