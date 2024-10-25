function data = getsinglecellmasks(data,reRunBgCell)
%Get all masks for single cell analysis;

%% try to get loaded information
try
    mask = data.mask;
catch
    mask = [];
end

bgmask = [];
if reRunBgCell==0
    try
        bgmask = data.bgmask;
    catch
        bgmask = [];
    end
end

try
    numCells = data.numCells;
catch
    numCells = 0;
end

%%open file
R2=bfopen(data.Location);
pics=R2{1};
pics=pics(:,1);
pn = length(pics);

for i=1:pn
    IMG(:,:,i)=pics{i};
end
images=double(IMG);
sx=size(images,1);
sy=size(images,2);
totalintensityimg = sum(images,3);

clear pics
clear IMG

%% try to get time data
try
    for i=1:pn
        T(i)=R{4}.getPlaneDeltaT(0, i-1).value;
    end
catch
    T=0:0.5:pn;
end
T = double(T);

%% MASKING IMAGES
if isempty(mask)
    if numCells == 0
        inputvalue = input('How many cells?');
    else
        inputvalue = numCells;
    end
    figure
    imshow(mat2gray(totalintensityimg));
    mask=zeros(sx,sy);
    for i=1:inputvalue
        individualmask=zeros(sx,sy);
        disp('select cell')
        bf=imfreehand();
        useMask=createMask(bf);
        mask=mask+useMask*i;
        %allmasks(:,:,i) = individualmask+useMask;
        Area_i(i)=size(nonzeros(useMask),1);
    end
else
    figure
    overlay = imoverlay(mat2gray(totalintensityimg), mask, 'colormap', 'jet', 'facealpha', .75, 'ZeroAlpha', 0);
    imshow(overlay)
    Area_i = data.Area_i; 
end
%% choose background cell

if isempty(bgmask)
    answer = questdlg('Choose a background cell','Notification','OK','Cancel');
    bgmask = zeros(sx, sy);
    figure
    imshow(mat2gray(totalintensityimg));
    disp('Choose backgroundcell');
    bf=imfreehand();
    useMask=createMask(bf);
    bgmask=bgmask+useMask;
    Area_bg=size(nonzeros(useMask),1);
else
    figure
    overlay = imoverlay(mat2gray(totalintensityimg), bgmask, 'colormap', 'jet', 'facealpha', .75, 'ZeroAlpha', 0);
    imshow(overlay)
    Area_bg = data.Area_bg; 
end

data.Area_i = Area_i;
data.mask = mask;
data.Area_bg = Area_bg;
data.bgmask = bgmask;
data.images = images;
data.T = T;

figure
imagesc(mask);
colorbar;

%data.allmasks=allmasks;
%data.activitycutoff=activitycutoff;

end

