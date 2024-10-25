function [NoSignalMask, NoSignalArea]=RemovingAreaswNoSignal(peakAmpMap,thresh)

[sx,sy]=size(peakAmpMap);
NoSignalMap=zeros(sx,sy);
NoSignalMask=zeros(sx,sy);
%%Nurin's original threshhold cutoff from averages of no signal area

NoSignalMap(find(peakAmpMap>=0.001 & peakAmpMap<=thresh*(1)))=1;
%  NoSignalMask=logical(NoSignalMap);
NoSignalMap=imfilter(NoSignalMap, [3 3], 'replicate');
NoSignalMap=bwlabeln(NoSignalMap,8);
NoSignalMap(find(peakAmpMap==0))=0;
NoSignalMask(find(NoSignalMap>0))=1;
NoSignalMask=logical(NoSignalMask);

NoSignalArea=length(nonzeros(NoSignalMap));
%   NoSignalArea=sum(stats)

  
   %% plot
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
t = linspace(0,1,n)';
r = interp1(t0,colors(:,1),t);
g = interp1(t0,colors(:,2),t);
b = interp1(t0,colors(:,3),t);
cmap = [r,g,b];

NoSignalMap(find(NoSignalMap==0))=NaN;

% imoverlay(mat2gray(imagesRaw(:,:,1)),NoSignalMap,[],[],cmap,0.3)

%imoverlayNurin(mat2gray(peakAmpMap),NoSignalMap,[],[],cmap,0.3);
%ButtonName = questdlg('Verify no signal areas look appropriate',...
 %                             'Alert','Ok','Cancel','Cancel');