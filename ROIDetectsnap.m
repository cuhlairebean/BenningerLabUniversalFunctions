%%ROIDETECT FOR SNAPS


SelImg = uigetfile('*.*');

R = bfopen(SelImg);
% ImAv = mean(R);
% HSV = ones(256,256,3);
% ImAvn = ImAv/max(ImAv(:));
% HSV(:,:,3) = ImAvn.^0.8;
% RGB = hsv2rgb(HSV);
% % RGB = single(RGB);
% figure
% imshow(RGB);

% Getting rid of Background
fig = figure('Name','Background');
ImGray = rgb2gray(R);
imagesc(ImGray)
disp('Select Background Signal')
bkgrndC=createMask(imfreehand);
close(fig);
BackC=ImGray.*bkgrndC;
BackC(~logical(BackC))=nan;
thresh = nanmax(nanmax(BackC));

[ImgNon, ~] = RemovingAreaswNoSignal(ImGray,thresh);
ImGray(ImgNon)=0;
figure
imshow(ImGray);

% SmRemv = RemovingSmallAreas(ImGray,10);
% figure
% imshow(SmRemv);
SmRemv = single(ImGray);
edgeThresh = 0.1;
amount = 1;
figure;
ImLoCo = localcontrast(SmRemv,edgeThresh,amount);
imshow(ImLoCo);

figure
ImShp = imsharpen(ImLoCo,'Radius',0.5,'Amount',5,'Threshold',0.1);
imshow(ImShp);

Thresh = 0.15;

BWCan = edge(ImShp,'Canny',Thresh);
% BWPrew = edge(IMG,'Prewitt');
% BWSob = edge(IMG,'Sobel');
% BWRob = edge(IMG,'Roberts');

% se90 = strel('line',2,90);
% se0 = strel('line',2,0);
% BWdil = imdilate(BWCan,[se0 se90]);
% se2 = strel('disk',1);
% BWop = imopen(BWdil,se2);

% BWfin = BWop;
% figure(5)
% imshow(BWfin);

% Overlay = imfuse(RGB,BWCan,'falsecolor','ColorChannels',[2 1 0]);
% figure(5);
% imagesc(Overlay);

BWoutline = bwperim(BWCan);
Segout = RGB;
Segout(BWoutline) = 255;
Final = imoverlay(RGB,BWoutline,'white');
figfin = figure;
imshow(Final);

output_dir = 'G:\BenningerLab\VKIsletROIOutput\';
prompt = {'Input File Name for Saving'};
title = 'Input';
dims = [1 35];
definput = {''};
file = inputdlg(prompt,title,dims,definput);
saveas(figfin,[output_dir file],'.tif');