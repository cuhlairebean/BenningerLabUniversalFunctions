function [ActiveArea, ActMap, ActMask, CoorMap] = getactivearea(firstimage, LogicalMap, peakAmpMap, coeffcorrMap, thresh)
%% CALCULATING ACTIVITY and GENERATING ACTIVITY MAPS
    [sx, sy] = size(peakAmpMap);
    ActMap=zeros(sx,sy);
    %Determine the active area using the bachground threshhold

    ActMap(peakAmpMap>=1.5*thresh)=1;

    
    ActMap(~LogicalMap)=0;
    ActMask=logical(ActMap);
    ActMap(~ActMap)=NaN;

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
    
    %%%%%%Use the CoorMap to overlay the activity map to put a number to activityremoving regions in the no signal mask, outside of the activity mask, and the image mask
    CoorMap=coeffcorrMap;
    CoorMap(~ActMask)=0;
    
    %CoorMap(CoorMap<.15)=0;  %%Another threshhold that is useful sometimes especially with movement
    ActiveArea=size(nonzeros(CoorMap),1);
    
        %% PLOTTING ACTIVITY MAP
    CoorMap(~logical(CoorMap))=NaN;
    pause(2)
    imoverlayNurin(mat2gray(firstimage(:,:,1)),CoorMap,[0, 1],[],cmap2,0.3);
    colorbar
    caxis([0 1])
    title('Activity Map')
end

