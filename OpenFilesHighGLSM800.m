function [ Dat ] = OpenFilesHighGLSM800(dir,files, mask)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

direc=1;
R=bfopen(dir);
if isempty(mask)
    mask=0;
end
%mask=files.mask;

% if length(find(mask))<1000
%     mask=0;
% end

pics=R{1};
pics=pics(:,1);
if isempty(R{4}.getPlaneDeltaT(0,1))
    T=1:size(pics,1);
else
    for i=1:length(pics)
        T(i)=R{4}.getPlaneDeltaT(0,i-1).value;
    end
end

T=double(T);
for i=1:length(pics)
    IMG(:,:,i)=pics{i};
end
pics={};

[sX,sY]=size(IMG(:,:,1));

stimMask=zeros(sX,sY);

ROICount=floor(double(R{4}.getROICount()/2));


% if TimeSplit==0
%     IMG=IMG(:,:,1:floor(size(IMG,3)/2));
%     T=T(1:floor(length(T)/2));
% else
%     IMG=IMG(:,:,floor(size(IMG,3)/2):end);
%     T=T(floor(length(T)/2):end);
% end
if isempty(files) % get stimulated area
    
    x=abs(floor(double(R{4}.getRectangleX(0,0))));
    y=abs(floor(double(R{4}.getRectangleY(0,0))));
    
    width=floor(double(R{4}.getRectangleWidth(0,0)));
    height=floor(double(R{4}.getRectangleHeight(0,0)));
    AVal=ones(height, width);
    stimMask(y:1:(y+height-1),x:1:(x+width-1))=AVal;
    
    %stimMask(sY-(height+1:1:y),x:1:(x+width-1))=AVal;
    %stimMask = flipud(stimMask);
    figure
    overlay = imoverlay(mat2gray(IMG(:,:,1)), stimMask, 'colormap', 'jet', 'facealpha', .75, 'ZeroAlpha', 0);
    imshow(overlay);
else
    
    stimMask=files.mask;
    
end
%[Dat]=HIanalyzeV3(IMG,T,0,direc);
[Dat]=HIanalyzeHighGLSM800(IMG,T,0,direc,stimMask,mask);

Dat.stimMask=stimMask;


Dat.img=IMG(:,:,1);
Dat.T=T;

%show correlation map
%% SET MAP COLORS
colors = [
    1 0 1 % First element = purple
    0 0 1 % blue
    0 1 0 % green
    1 1 0 % yellow
    1 .65 0 % orange
    1 0 0]; % red

n = 256; % size of new color map
m = size(colors,1);
t0 = linspace(0,1,m)';
tt = linspace(0,1,n)';
r = interp1(t0,colors(:,1),tt);
g = interp1(t0,colors(:,2),tt);
b = interp1(t0,colors(:,3),tt);
cmap2 = [r,g,b];

newcorrMap = Dat.corrMap;
% %newcorrMap(find(isnan(newcorrMap)))=0;
% 
% overlay = imoverlay(mat2gray(IMG(:,:,1)), newcorrMap, 'facealpha', .75, 'ZeroAlpha', 0);
% imshow(overlay);
pause(2)
imoverlayNurin(mat2gray(IMG(:,:,1)), newcorrMap,[0, 1],[],cmap2,0.3);
colorbar
caxis([0 1])

%save as tiff
imagewd = getframe(gcf);
savedfile = regexprep(dir,'.czi','.tiff');
imwrite(imagewd.cdata, savedfile, 'Resolution',300);
end

