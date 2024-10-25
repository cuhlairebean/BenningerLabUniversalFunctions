function [ data ] = GetMask_JD(data)


images=double(data.images);
imagesRaw=mat2gray(images);
%%ImgPrep
if data.filter ~=1
    matrixClusterLo=imagesRaw;
    h = fspecial('average', 10);
    images=imfilter(matrixClusterLo,h,'replicate');
else
    images = imagesRaw;
end

mask = data.mask;
timeVec = data.TSize;

%xrealigns images if Islets are moving
%images=registerImages(images);
imageMean=mean(images,3);
img=images(:,:,1);
sx=size(img,1);
sy=size(img,2);
sz=size(images,3);
L=sz;
NFFT=2^nextpow2(L);
Fs=mean(diff(timeVec));
Fs=1/Fs;
TF=Fs./2*linspace(0,1,NFFT/2+1);
TF(1)=[];
TF=1./TF;

%Mask first image
if isempty(mask);
    imshow(mat2gray(imagesRaw(:,:,6)));
    %imagesc(matrixMask)
    %numAreas=input('Enter number of areas');
    numAreas=1;
    mask=zeros(sx,sy);

    % level=graythresh(mat2gray(images(:,:,1)));
    % mask=im2bw(mat2gray(images(:,:,1)),level);
    for i=1:numAreas
        matrixDrawMask = imfreehand;
        ButtonName = questdlg('Do you want to proceed?',...
                              'Question');
        if strcmp(ButtonName, 'Yes')
            useMask = createMask(matrixDrawMask);
        else
            %break
        end
        data.mask=mask+useMask;
    end
end
%labelbw=bwlabel(mask,4);
%numAreas=max(labelbw(:));
%currArea=1;
for i=1:sz
    if data.filter ==1
        %filters noisier data...
        images(:,:,i) = medfilt2(images(:,:,i), [10,10]).*data.mask;
    else 
        images(:,:,i) = images(:,:,i).*data.mask;
    end
end 
data.images = images;

end
