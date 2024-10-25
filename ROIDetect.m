close all
clear all
clc
%% Loading the Ca Imaging CZI
output_dir = 'G:\BenningerLab\VKIsletROIOutput\';
SelImg = uigetfile('*.*','Select Image');

R = bfopen(SelImg);

if length(R(:,4))>1
    pics=R(:,1);
    for i=1:length(pics)
        pics1 = pics{i};
        pics2{i} = pics1{1};
    end
    pics=pics2';
else
    pics=R{1};
    pics=pics(:,1);
end

% try R{4}.getPlaneDeltaT(0,2)
%     
%     if ~isempty(R{4}.getPlaneDeltaT(0,0))
%         for i=1:length(pics)
%             T(i)=R{4}.getPl/aneDeltaT(0,i-1);
%         end
%     else
%         T=1:size(pics,3);
%     end
%     
% catch
%     
%     % for i=1:length(pics)
%     % T(i)=R{4}.getPlaneDeltaT(i-1,0);
%     % end
%     
%     T=0:1:length(pics);
% end

% T=double(T);

for i=1:length(pics)
    IMG(:,:,i)=pics{i};
end
pics={};

Images=double(IMG);
RawImg=Images(:,:,1);


%% CODE TO MODIFY TCs
st=400;
%ed=size(Images,3);
ed=800;
T=st:1:ed;
%
images=Images(:,:,st:ed);
DataOut.images=images;
% % % %

%% DECLARING IMAGE PROPERTIES

sx=size(images,1);
sy=size(images,2);
sz=size(images,3);


for i=1:sz
    %images1(:,:,i)=imfilter(imagesRaw(:,:,i),h,'replicate').*mask;
    %images2(:,:,i)=imagesRaw(:,:,i).*mask;
    images(:,:,i)=medfilt2(images(:,:,i),[5 5]);
end

% series1 = data{1, 1};
% mov1 = cat(series1(1,:));
% cmap = [zeros(256,1), linspace(0,1,256)', zeros(256,1)]; % Green Colormap
% MOV = immovie(series1{1:1426,:},cmap);
% % colormap(cmap);
% % colorbar;
% figure(1)
% implay(MOV);

ImAv = mean(images,3);
% figure
% imagesc(ImAv);
HSV = ones(512,512,3);
HSV(:,:,1) = 0.3;
ImAvn = ImAv/max(ImAv(:));
HSV(:,:,3) = ImAvn.^0.8;
RGB = hsv2rgb(HSV);
% RGB = single(RGB);
figure
imshow(RGB);

% Getting rid of Background
fig = figure('Name','Background');
ImGray = rgb2gray(RGB);
imagesc(ImGray)
disp('Select Background Signal')
bkgrndC=createMask(imfreehand);
close(fig);
BackC=RGB.*bkgrndC;
BackC(~logical(BackC))=nan;
thresh = nanmax(nanmax(BackC));

[ImgNon, ~] = RemovingAreaswNoSignal(RGB,thresh);
RGB(ImgNon)=0;
% figure
% imshow(RGB);

% SmRemv = RemovingSmallAreas(ImGray,10);
% figure
% imshow(SmRemv);
SmRemv = single(RGB);
edgeThresh = 0.1;
amount = 1;
ImLoCo = localcontrast(SmRemv,edgeThresh,amount);
% figure;
% imshow(ImLoCo);

FirstFig = figure('Name','Img without Outline');
ImShp = imsharpen(ImLoCo,'Radius',0.5,'Amount',5,'Threshold',0.1);
imshow(ImShp);

% cmap = [zeros(256,1), linspace(0,1,256)', zeros(256,1)]; % Green Colormap
% ImShp = uint16(ImShp);

Thresh = 0.15;

% Just fucking around with binarizing for singling out disconnected bits
% % ImGray = rgb2gray(ImShp);
% % ImBW = imbinarize(ImGray);
% % figure;
% % imshow(ImBW);
% % 
% % SE1 = strel('disk',4);
% % ErodIm = imerode(ImBW,SE1);
% % figure;
% % imshow(ErodIm);
% % 
% % SE2 = strel('diamond',1);
% % DilErodIm = imdilate(ErodIm,SE2);
% % figure;
% % imshow(DilErodIm);

BWCan = edge(ImGray,'Canny',Thresh);
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

% Applies edge detection to give rough outline of cell boundaries 
BWoutline = bwperim(BWCan);
Segout = ImShp;
Segout(BWoutline) = 255;
Final = imoverlay(ImShp,BWoutline,'k');
FigFinal = figure('Name','Img with Outline');
imshow(Final);

% For Selecting multiple ROIs
n = 1;
j = 1;

% while n == 1
%     ROImask = drawfreehand();
%     ROImask = createMask(ROImask);
%     ROI(j) = bwboundaries(ROImask);
%     ROIx{j} = ROI{j,j}(:,2);
%     ROIy{j} = -1*ROI{j,j}(:,1);
%     set(gca,'XTick',[],'YTick',[]); %Removes axis tick marks from image
%     UImp = input('Select Additonal ROI? 1=Yes, 0=No \n');
%     if UImp == 1
%         n = 1;
%         j = j+1;
%     elseif UImp == 0
%         n = 0;
%     end
% end

while n == 1
    ROImask = drawfreehand();
    ROImask = createMask(ROImask);
    ROImask(:,:,j) = ROImask; % Need to figure out why it keeps overwriting the first page of the 3D array...
    set(gca,'XTick',[],'YTick',[]); %Removes axis tick marks from image
    UImp = input('Select Additonal ROI? 1=Yes, 0=No \n');
    if UImp == 1
        n = 1;
        j = j+1;
    elseif UImp == 0
        n = 0;
    end
end
    



prompt = {'Input File Name for Saving'};
title = 'Input';
dims = [1 35];
definput = {''};
file = inputdlg(prompt,title,dims,definput);
filename = [output_dir file{1}];
saveas(FigFinal,filename,'tiffn');
saveas(FirstFig,[filename '_no_outline'],'tiffn');

close all
msgbox('Images Saved');
