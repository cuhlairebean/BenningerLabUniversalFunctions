function [DataOut]=GFPinCaps(data)
%% LOADING IMAGES
DataOut = data;
if length(fieldnames(data))==1
    backgroundmask=[];
else
    backgroundmask = data.backgroundmask;
end
R=bfopen(data.Location);
pics=R{1};
pics=pics(:,1);
pn = length(pics);

for i=1:pn
    IMG(:,:,i)=pics{i};
end

clear pics

Images=double(IMG);
st=1;
ed=1;

images=Images(:,:,st:ed);
DataOut.images=images;
Images=double(images);
RawImg=Images(:,:,1);
OrigImag = RawImg;
imagesRaw=images;

sx=size(images,1);
sy=size(images,2);
sz=ed;

%h = fspecial('average', 5);
for i=1:sz
    %images1(:,:,i)=imfilter(imagesRaw(:,:,i),h,'replicate').*mask;
    %images2(:,:,i)=imagesRaw(:,:,i).*mask;
    images(:,:,i)=medfilt2(imagesRaw(:,:,i),[5 5]);
end
% 
% BW = imbinarize(mat2gray(images(:,:,1)),'adaptive','Sensitivity',0.6);
% figure;
% imshow(BW)
% BW=medfilt2(BW,[5 5]);
% 
% BW2 = imfill(BW,'holes');
% figure
% imshow(BW2)
% title('Filled Image')

%% GFP image as mask %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(backgroundmask)
    figure;
    imagesc(RawImg)
    disp('select background area')
    backgroundmask=createMask(imfreehand);
else
    figure
    overlay = imoverlay(mat2gray(RawImg),backgroundmask, 'colormap', 'jet', 'facealpha', .75, 'ZeroAlpha', 0);
    imshow(overlay)
end
backgroundMap=RawImg.*backgroundmask;
thresh=max(max(backgroundMap));
[nosignalmask, ~]=RemovingAreaswNoSignal(RawImg,thresh);
RawImg(nosignalmask)=0;

figure
overlay = imoverlay(mat2gray(OrigImag), RawImg, 'colormap', 'jet', 'facealpha', .75, 'ZeroAlpha', 0);
imshow(overlay)
colorbar


BW = imbinarize(mat2gray(images(:,:,1)),'adaptive','Sensitivity',0.6);
figure;
imshow(BW)
BW=medfilt2(BW,[5 5]);

BW2 = imfill(BW,'holes');
figure
imshow(BW2)
title('Filled Image')

intmap=bwlabeln(im2bw(RawImg));

figure
imagesc(im2bw(intmap))
title('Map of Areas')

%s = regionprops(BW, I, {'Centroid','WeightedCentroid'});

STATS=regionprops(intmap,'area');
DataOut.STATS=STATS;
for j=1:length(STATS)
    stats(j)=STATS(j).Area;
end

GFPCapsMap=zeros(size(map));
nostats=find(stats<500);
for i=1:length(nostats)
    GFPCapsMap(intmap==nostats(i))=1;
end
GFPCapsMap = logical(GFPCapsMap);
pause(2)

CapImage = RawImg;
NotCapsImage = RawImg
CapImage(~GFPCapsMap)=0;
NotCapsImage(GFPCapsMap)=0;

figure
overlay = imoverlay(mat2gray(OrigImg), CapsImage, 'colormap', 'jet', 'facealpha', .75, 'ZeroAlpha', 0);
imshow(overlay)
colorbar

DataOut.backgroundmask=backgroundmask;
DataOut.thresh=thresh;
end
