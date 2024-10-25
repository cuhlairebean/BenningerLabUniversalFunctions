function [NoSignalMask, NoSignalArea]=RemovingAreaswNoSignal2(avgintmap,thresh,peakAmpMap,thresh2)

[sx,sy]=size(avgintmap);

NoSignalMap=zeros(sx,sy);
NoSignalMask=zeros(sx,sy);

NoSignalMap2=zeros(sx,sy);

   % NoSignalMap(find(avgintmap>=0.001 & avgintmap<=thresh*(1+0.80)))=1;
   % %Nurin's original threshhold cutoff from averages of no signal area
   NoSignalMap(find(avgintmap<=thresh & peakAmpMap<=2*thresh2))=1;
  % NoSignalMap2(find(peakAmpMap>thresh2*(2)))=1;
  % NoSignalMap(find(NoSignalMap2==1))=0;
   %  NoSignalMask=logical(NoSignalMap);
   NoSignalMap=imfilter(NoSignalMap, [2 2], 'replicate');
   NoSignalMap=bwlabeln(NoSignalMap,8);
   NoSignalMap(find(avgintmap==0))=1;
   NoSignalMask(find(NoSignalMap>0))=1;
   NoSignalMask=logical(NoSignalMask);
%%
STATS=regionprops(NoSignalMap,'area');


%    for j=1:length(STATS)
%       stats(j)=STATS(j).Area; 
%    end
%    
%    nostats=find(stats<20);
%    yesstats=find(stats>5);
   %%
%    NoSignalArea=sum(stats)
   
% for i=1:length(nostats)
% NoSignalMap(find(NoSignalMap==nostats(i)))=0;
% end
% 
% for i=1:length(yesstats)
%     NoSignalMap(find(NoSignalMap==yesstats(i)))=1;
% end
% % 
% a=sum(yesstats)
 
  NoSignalArea=length(nonzeros(NoSignalMap));
%   NoSignalArea=sum(stats)
  
%   
%    %% plot
%      colors = [
% 1 0 1 % First element = purple
% 0 0 1 % blue
% 0 1 0 % green
% 1 1 0 % yellow
% 1 .65 0 % orange
% 1 0 0]; % red
% 
% n = 256; % size of new color map
% m = size(colors,1);
% t0 = linspace(0,1,m)';
% t = linspace(0,1,n)';
% r = interp1(t0,colors(:,1),t);
% g = interp1(t0,colors(:,2),t);
% b = interp1(t0,colors(:,3),t);
% cmap = [r,g,b];
% 
% 
%    NoSignalMap(find(NoSignalMap==0))=NaN;
% 
% %    imoverlay(mat2gray(imagesRaw(:,:,1)),NoSignalMap,[],[],cmap,0.3)
% %
% imoverlayNurin(mat2gray(avgintmap),NoSignalMap,[],[],cmap,0.3);
