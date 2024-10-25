%% ROI TRACING PROGRAM
% Credit: W. Schleicher, April 2020
% Designed to import a calcium trace file for which ROIs have already been
% drawn (in FIJI) and import a representative image on which the ROIs and
% labels are visible. The user then traces these ROIs to be saved and later
% imported into another program.
% NOTE: Doesn't work on calcium trace videos that include multiple islets,
% unless the representative image with ROIs is of the same magnification.

% NOTE: When tracing ROIs, no mistakes can be made!!! If you mess up, you
% have to start over. Otherwise the mask set will not import to the
% analysis program correctly.

close all
clear all
clc

addpath('H:\SCHLEICHER\MatLabCode\UniCode\');
%%

R = bfopen('Select Calcium Imaging File');

pics=R{1};
pics=pics(:,1);

for i=1:length(pics)
    IMG(:,:,i)=pics{i};
end

images=double(IMG); % converts images to double precision
RawImg=images(:,:,1); % assigns the first frame of the video to RawImg variable

clear pics R IMG;
clear omeMeta;

[xx, yy, zz] = size(images(:,:,:));

[file,path] = uigetfile('*.*','Select Image with ROIs');
addpath(path);
OGFig = imread(file);
OGFig = imresize(OGFig,[xx yy]);

disp('Select Folder for Output Directory');
output_dir = uigetdir;
filename = input('Input Filename\n','s');

MaskFig = figure('Name','Trace around ROIs');
imshow(OGFig);

CellMask = double(zeros(xx,yy));


numcells = 1;
k = 1;
while k > 0
    disp('Trace Around ROIs')
    ROIMask = imfreehand(); %User draws region around cell
    ROIMask = createMask(ROIMask); %Mask is created from drawn region
    CellMask = CellMask + ROIMask.*numcells; %CellMask array is updated with new mask; new mask is multiplied by the cell label before updating
    CellMask(find(CellMask>numcells)) = numcells; %If a region is overlapped, it is instead attributed to the most recent region
    UImp = input('Select Additonal ROI? 1=Yes, 0=No \n'); % No preallocated number of ROIs, User input required to continue ROI drawing or to stop
    if UImp == 1 %User input to determine if another region is drawn
        k = k+1;
        numcells = numcells+1;
    elseif UImp ~= 1
        k = 0;
    end
end

save([output_dir '\' filename 'Masks' '.mat'],'CellMask'); %Saves CellMask array
save([output_dir '\' filename 'CellNumber' '.mat'],'numcells'); %Saves numcells array

disp('End');
