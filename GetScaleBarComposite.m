function []=GetScaleBarComposite(LocationGFP, LocationCa, savefilename,samefile, cachannel, imgnum, x_location, y_location, Scalebarlength) % in pixels

R=bfopen(LocationCa);
pics=R{1};
pics=pics(:,1);
pn = length(pics);

for i=1:pn
    IMG(:,:,i)=pics{i};
end
Images=double(IMG);
ImagesCa = Images;
clear pics IMG

if samefile
    ImagesCa = Images(:,:,cachannel:2:end);
    if cachannel ==1
        ImagesGFP = Images(:,:,2:2:end);
    else
        ImagesGFP = Images(:,:,1:2:end);
    end
else
    RGFP=bfopen(LocationGFP);
    pics=RGFP{1};
    pics=pics(:,1);
    pn = length(pics);
    
    for i=1:pn
        IMG(:,:,i)=pics{i};
    end
    ImagesGFP=double(IMG);
    clear pics IMG
end

RawImgCa=ImagesCa(:,:,imgnum);
RawImgGFP=ImagesGFP(:,:,imgnum);

HSV = ones(size(ImagesGFP,1),size(ImagesGFP,2),3);
HSV(:,:,1) = HSV(:,:,1)*.333;
%HSV for green is (120 degrees/360, 100%, 50%)
ImgInt = RawImgGFP/max(RawImgGFP(:));
HSV(:,:,3) = ImgInt;%.*1.5;
RGB = hsv2rgb(HSV);
RGB = single(RGB);
% figure
% h = imshow(RGB);
% hold on

HSV2 = ones(size(ImagesCa,1),size(ImagesCa,2),3);
HSV2(:,:,1) = HSV2(:,:,1)*25/360;
%HSV red(0°,100%,100%)
%HSV for green is (120 degrees/360, 100%, 50%)
ImgInt2 = RawImgCa/max(RawImgCa(:));
HSV2(:,:,3) = ImgInt2;%.*1.5;
RGB2 = hsv2rgb(HSV2);
RGB2 = single(RGB2);
% figure
% h = imshow(RGB2);


% C = imfuse(RGB, RGB2);
% C = imfuse(RGB, RGB2,'falsecolor','Scaling','joint','ColorChannels',[1 2 0]);
% figure
% C = imfuse(RGB, RGB2,'blend','Scaling','independent');
% imshow(C);

figure
C = imadd(RGB, RGB2);
imshow(C);

sx=size(ImagesCa,1);
sy=size(ImagesCa,2);
sz=size(ImagesCa,3);

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
C = quiver(x_location,y_location,Scalebar_pixels,0,'ShowArrowHead','off','Color','w','LineWidth',3.0);
txt = [num2str(Scalebarlength) '\mum'];
%text(x_location,y_location-20,txt,'Color','w');

imagewd = getframe(gcf);
imwrite(imagewd.cdata, [savefilename '.tiff'], 'Resolution',600);

end