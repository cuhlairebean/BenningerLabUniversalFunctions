%mask =[];    
bkgcell = [];

 %% MASKING IMAGES
    if isempty(mask)
        figure
        %imshow(mat2gray(mean(images,3)));
        imagesc(mean(images,3));
        %imshow(imagesc(mean(channel11,3)+mean(channel22,3)+mean(channel33,3)));
        numAreas=1;
        mask=zeros(sx,sy);
        title('Average Intensity')
        disp('select islet')
        for i=1:numAreas
            bf=imfreehand();
            useMask=createMask(bf);
            mask=mask+useMask;
        end
    else
        figure
        overlay = imoverlay(mat2gray(RawImg), mask, 'colormap', 'jet', 'facealpha', .75, 'ZeroAlpha', 0);
        imshow(overlay);
    end
    DataOut.mask(:,:,zz) = mask;
    imagesRaw = images;
    imagesfilt = smoothimage(imagesRaw, 5); %smooths images with 5X5
    
    images=imagesfilt.*mask;
    
    %% Peak Amp Map
    [peakAmpMap, ~, ~, coeffcorrMap] = getpeakampmap(T, images, mask);
    DataOut.PeakAmpMap(:,:,zz)=peakAmpMap;
    %    DataOut.StdDevMap = stddevMap;
    
    %% APPLYING "background" THRESHOLD
    AverageIntensityMap = mean(images, 3).*mask;
   % figure; imagesc(AverageIntensityMap);
    colorbar;
    if isempty(bkgcell)
        figure;
        imagesc(peakAmpMap)
        title('peakAmpMap')
        disp('select background cell')
        bkgcell=createMask(imfreehand);
    else
        figure
        overlay = imoverlay(mat2gray(RawImg), bkgcell, 'colormap', 'jet', 'facealpha', .75, 'ZeroAlpha', 0);
        imshow(overlay)
    end
    BKGMap=AverageIntensityMap.*bkgcell;
    BKGMap(~logical(BKGMap))=nan;
    bkgthresh=nanmax(nanmax(BKGMap));
    
    DataOut.BKGMap(:,:,zz)=BKGMap;
    DataOut.BKGthresh(zz)=bkgthresh;
    DataOut.bkgcell(:,:,zz)=bkgcell;
    
    %%%%%%%%%%%%silentcell%%%%%%%%%%%%%%%%%
    %%%% using the silent/background cell threshold to remove areas with NO signal %%%
    SCMap=peakAmpMap.*bkgcell;
    SCMap(~logical(SCMap))=nan;
    scthresh = nanmean(nanmean(SCMap));
    [NoSignalMask, ~]=RemovingAreaswNoSignal2(AverageIntensityMap,scthresh, DataOut.PeakAmpMap(:,:,zz), scthresh);
    
    SCMap=peakAmpMap.*NoSignalMask;
    SCMap(~logical(SCMap))=nan;
    nanmean(nanmean(SCMap))
    scthresh = nanmedian(nanmedian(SCMap));
    
    DataOut.scthresh(zz)=scthresh;
    
    figure
    overlay = imoverlay(mat2gray(RawImg), NoSignalMask, 'colormap', 'jet', 'facealpha', .75, 'ZeroAlpha', 0);
    imshow(overlay)
    title('Background Removed');
    
    %% %%%%%% Final Mask to use%%%%%%%%%%%%%%%%%%%%%%%%
    LogicalMap=zeros(sx,sy);
    LogicalMap(logical(mask))=1;
    LogicalMap(NoSignalMask)=0;
    TotalArea=length(nonzeros(LogicalMap));
    DataOut.TotalArea(zz) = TotalArea;
    Area=length(nonzeros(LogicalMap));
    LogicalMap = logical(LogicalMap);
    
    figure
    overlay = imoverlay(mat2gray(RawImg), LogicalMap, 'colormap', 'jet', 'facealpha', .75, 'ZeroAlpha', 0);
    imshow(overlay)
     title('All Area Analyzed');
    
    %% CALCULATING ACTIVITY and GENERATING ACTIVITY MAPS
    [ActiveArea, ActMap, ActMask, CoorMap] = getactivearea(RawImg, LogicalMap, peakAmpMap, coeffcorrMap, scthresh);
    
    DataOut.RatioActive(zz)=ActiveArea/(Area);%% CALCULATING PERCENT ACTIVE AREA
    DataOut.ActiveArea(zz) = ActiveArea;
    DataOut.Area(zz) = Area;
    disp(['Active Area = ' num2str(DataOut.RatioActive(zz))]);%print active are
    AreaStacks = AreaStacks + Area;
    ActiveAreaStacks = ActiveAreaStacks + ActiveArea;
    
    %% SAVING ACTIVITY VARIABLES
    %set a breakpoint and look at graphs
    DataOut.ActMask(:,:,zz)=ActMask;
    DataOut.ActMap(:,:,zz)=ActMap;
    DataOut.ActMapValues(:,:,zz)=CoorMap;
    DataOut.CorrMap(:,:,zz)=coeffcorrMap;
    DataOut.Location=data.Location;
    DataOut.NoSignalMask(:,:,zz)=NoSignalMask;
    DataOut.LogicalMap(:,:,zz)=LogicalMap;
    DataOut.mask(:,:,zz)=mask;
    DataOut.imagesfilt(:,:,:,zz)=images;