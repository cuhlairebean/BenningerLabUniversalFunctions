function []=GetScaleBar(Location,savefilename,samefile, whichchannel, imgnum, x_location, y_location, Scalebarlength, graygreenred) % in pixels

R=bfopen(Location);
pics=R{1};
pics=pics(:,1);
pn = length(pics);

for i=1:pn
    IMG(:,:,i)=pics{i};
end
clear pics
Images=double(IMG);

if samefile
    Images = Images(:,:,whichchannel:2:end);
end

%Can adjust the time frame of images here%
st=1;
ed=size(Images,3);
T=st:.5:ed;
t=st:ed; %for duty cycle

images=Images(:,:,st:ed);
RawImg=images(:,:,imgnum);

figure
if graygreenred==1
    h = imshow(mat2gray(RawImg));
elseif graygreenred==2
    %%green
    HSV = ones(size(Images,1),size(Images,2),3);
    HSV(:,:,1) = HSV(:,:,1)*.333;
    
    %HSV for green is (120 degrees/360, 100%, 50%)
    ImgInt = RawImg/max(RawImg(:));
    HSV(:,:,3) = ImgInt;%.*1.5;
    RGB = hsv2rgb(HSV);
    RGB = single(RGB);
    h = imshow(RGB);
elseif graygreenred==3
    %%red
    HSV = ones(size(Images,1),size(Images,2),3);
    HSV(:,:,1) = HSV(:,:,1)*25/360;
    %HSV red(0°,100%,100%)
    %HSV for green is (120 degrees/360, 100%, 50%)
    ImgInt = RawImg/max(RawImg(:));
    HSV(:,:,3) = ImgInt;%.*1.5;
    RGB = hsv2rgb(HSV);
    RGB = single(RGB);
    h = imshow(RGB);
    %imshow(mat2gray(RawImg));
    %cmap = [zeros(256,1), linspace(0,1,256)', zeros(256,1)]; % Green Colormap
    %imagesc(RawImg);
    %colormap(cmap);
end


sx=size(images,1);
sy=size(images,2);
sz=size(images,3);

%%get physical pixel size

omeMeta = R{1, 4};
stackSizeX = omeMeta.getPixelsSizeX(0).getValue(); % image width, in pixels
stackSizeY = omeMeta.getPixelsSizeY(0).getValue(); % image height, in pixels

voxelSizeXdefaultValue = omeMeta.getPixelsPhysicalSizeX(0).value(); % returns value in default unit
voxelSizeXdefaultUnit = omeMeta.getPixelsPhysicalSizeX(0).unit().getSymbol(); % returns the default unit type
voxelSizeX = omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER); % in µm
voxelSizeXdouble = voxelSizeX.doubleValue(); % The numeric value represented by this object after conversion to type double
voxelSizeY = omeMeta.getPixelsPhysicalSizeY(0).value(ome.units.UNITS.MICROMETER); % in µm
voxelSizeYdouble = voxelSizeY.doubleValue(); % The numeric value represented by this object after conversion to type double

PixelSize = voxelSizeXdouble; % in um
Scalebar_pixels = Scalebarlength/PixelSize; % in pixels

hold on
h = quiver(x_location,y_location,Scalebar_pixels,0,'ShowArrowHead','off','Color','w','LineWidth',3.0);
txt = [num2str(Scalebarlength) '\mum'];
%text(x_location,y_location-20,txt,'Color','w');

imagewd = getframe(gcf); 
imwrite(imagewd.cdata, [savefilename '.tiff'], 'Resolution',600);

end