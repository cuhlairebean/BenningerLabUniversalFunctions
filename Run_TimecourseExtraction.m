function [SingleCellData] = Run_TimecourseExtraction (filename, st, ed)

R = bfopen(filename); % Uses bfopen program to open .czi/.lsm image files

pics=R{1};
pics=pics(:,1);

for i=1:length(pics)
    IMG(:,:,i)=pics{i};
end

images=double(IMG); % converts images to double precision

images = images(:,:,1:1:end); % assuming only 1 channel 
RawImg=images(:,:,1); % assigns the first frame of the video to RawImg variable
ref = images(:,:,1:end); % This chooses from a set of frames to create the max intensity map: here 1:end
MaxIntensityMap = max(ref, [],3);
ref_pic = (mat2gray(MaxIntensityMap));

clear pics R IMG;
clear omeMeta;

%% CODE TO MODIFY TIME COURSES
% This block converts the 3rd dimension of array into timecourse

%st=1; % STARTING FRAME: Change this per sample if needed

[xx, yy, zz] = size(images(:,:,:));

%ed=zz; % ENDING FRAME: Change this per sample if needed

T=st:1:ed;

images=images(:,:,st:ed);

% DECLARING IMAGE PROPERTIES
sx=size(images,1);
sy=size(images,2);
sz=size(images,3);

for i=1:sz
    images(:,:,i)=medfilt2(images(:,:,i),[5 5]); %applies filter to clean up images
end

ImAv = mean(images,3); %compresses all frames into single array of intensities
HSV = ones(sx,sy,3); %preallocates a 3 dimensional array
ImAvn = ImAv/max(ImAv(:));
HSV(:,:,3) = ImAvn.^0.8; %evens out intensity across images
HSV(:,:,1) = 0.3333;%converts image to green image
RGB = hsv2rgb(HSV); %converts to rgb image: this might not be useful

OGFig = figure(1);
imshow(RGB);

%% 1. Single Cell Masking
% User Draws ROIs around each cell within the islet
% ROIs are saved in "CellMask" array and called back to throughout analysis
% Follow the prompts in the Command Window

% Getting rid of Background using JD code
GrayFig = figure('Name','Remove Background');
ImGray = rgb2gray(RGB); %converts image to gray
imagesc(ImGray) %displays in default colormap; easier to see cell borders
disp('Select Background Signal')
bkgrndC=createMask(imfreehand); %allows user to select background area to be removed
close(GrayFig);
BackC=ImGray.*bkgrndC; %multiplies mask by the image
BackC(~logical(BackC))= nan; %takes inverse of masked image and sets it to nan
thresh = nanmax(nanmax(BackC)); %sets threshold based on the max of values after nan's are removed
close(OGFig);

[ImgNon, ~] = RemovingAreaswNoSignal(ImGray,thresh); %Uses function from JD to remove background based on threshold
ImGray(ImgNon)=0; %sets the removed areas to 0
NoSigFig = figure('Name','Background Removed; Draw ROI around all cells');
imagesc(ImGray);
clear OGFig GrayFig;

% MASKING ISLET
% User draws ROI around islet to plot whole islet timecourse
disp('Draw ROI around entire area of interest');
h = imrect(gca); % try using a rectangle since we are looking at tons of cells
        pos_rect = h.getPosition();
        pos_rect = round(pos_rect);
        ROIMask_Start = images(pos_rect(2) + (0:pos_rect(4)), pos_rect(1) + (0:pos_rect(3)), :);
%ROIMask_Start = createMask(ROIMask_Start);
sz_FR=ed-st+1;


StartMaskStack = ROIMask_Start;
%StartMaskStack = images.*ROIMask_Start;
IsletTC = shiftdim(mean(mean(StartMaskStack)),2);
IsletTCfig = figure('Name','Whole Islet Time Course');
plot(IsletTC); % see if there is any general pattern 

CellMask = double(zeros(xx,yy));
numcells = 1;

answer = questdlg('Import Mask Set or Draw New Masks?', ...
    'Import Dialogue', ...
    'Import Mask Set','Draw New Masks','Draw New Masks');

switch answer
    case 'Import Mask Set'
        disp('Importing Mask Set')
        disp('Please Select Mask File');
        [Mask,mpath] = uigetfile('*.mat','Please select Mask file'); % User selects .mat mask file - only valid for MatLab generated files!!!
        addpath(mpath);
        dummyMask = struct2cell(load(Mask)); % Loads in Mask file
        CellMask = dummyMask{:,:}; % Assigns loaded mask file to CellMask Variable
        disp('Please Select Cell Number File');
        [num,npath] = uigetfile('*.mat','Please select Cell Number file',mpath); % User selects .mat cell number file - only valid for MatLab generated files!!!
        addpath(npath);
        dummyNum = struct2cell(load(num)); % Loads in Cell Number file
        numcells = dummyNum{:,:};
    case 'Draw New Masks'
        disp(answer)
        NoSigFig = figure('Name','Draw ROIs Around Cells of Interest');
        imshow(ref_pic);
        
        k = 1;
        while k > 0
            disp('Draw ROIs Around Cells of Interest')
            ROIMask = imfreehand(); %User draws region around cell
            ROIMask = createMask(ROIMask); %Mask is created from drawn region
            CellMask = CellMask + ROIMask.*numcells; %CellMask array is updated with new mask; new mask is multiplied by the cell label before updating
            CellMask(find(CellMask>numcells)) = numcells; %If a region is overlapped, it is instead attributed to the most recent region
            UInp = input('Select Additonal ROI? 1=Yes, 0=No \n'); % No preallocated number of ROIs, User input required to continue ROI drawing or to stop
            if UInp == 1 %User input to determine if another region is drawn
                k = k+1;
                numcells = numcells+1;
            elseif UInp ~= 1
                k = 0;
            end
        end
        
        save([output_dir '/' sample '_' 'Masks' '.mat'],'CellMask'); %Saves CellMask array
        close(NoSigFig);
        clear NoSigFig;
        
        save([output_dir '/' sample '_' 'CellNumber' '.mat'],'numcells'); %Saves numcells array
end

% TIMECOURSE
sz_2=ed-st; % includes all frames for exporting of entire trace

PlotHandles = zeros(1,numcells);
PlotLabels = cell(1,numcells);
CellTC = zeros(sz_2,numcells);
TC = zeros(sz_2,1);

tic
for i = 1:numcells
    TCMask = CellMask; %Pulls in CellMask array
    TCMask(find(TCMask ~= i)) = 0; %Gets rid of all masks besides current one
    MaskedIMGstack = images.*logical(TCMask); %Applies current mask to all frames of timecourse
    %get rid of for loop
    for ii = 1:sz_2 
        TCnoZero = MaskedIMGstack(:,:,ii); %Pulls in current frame
        TCnoZero = TCnoZero(TCnoZero>0); %Accounts for any zeros from preallocation
        TC(ii) = mean(TCnoZero); %Calculates mean intensity of current frame
    end

    CellTC(:,i) = TC; %Updates CellTC array with mean intensities from each frame of each cell
    PlotLabels{i} = [' Cell ' num2str(i)]; %Updates labels with current cell number for legend
end
toc

% Plotting traces for entire time course
TCFig = figure('Name','Average Intensity Over Time');
plot(CellTC);
legend(PlotLabels);
clear MaskedIMGstack;
saveas(TCFig,[output_dir '/' sample '_' 'TCPlot' '.tif']); %Saves figure of each cell's timecourse
save([output_dir '/CaWaveForm.mat'],'CellTC'); % same cawaveform - Cell TC

%% View Cells w labels
ImgWave = rgb2gray(RGB); %Calls in representative image in grayscale
Cent = zeros(numcells,2);
for j = 1:numcells
    DummyBW = CellMask; %Pulls in CellMask Array
    DummyBW(find(DummyBW ~= j)) = 0; %Gets rid of all but current mask
    STATS = regionprops(logical(DummyBW),'Centroid'); %Finds centroid of mask
    Cent(j,:) = STATS.Centroid; %Updates array with current centroid coordinates
end

% show figure
cell_perims = bwperim(CellMask);
cells_outline = imfuse(ref_pic, cell_perims);
cells_w_labels = figure; 
imshow(cells_outline)
title('Whole Islet')
ci = 1;
for c = 1:numcells
    text(Cent(ci,1),Cent(ci,2),num2str(c),'Color','w'); %Labels cells in the image with their respective region number
    ci = ci+1;
end


saveas(cells_w_labels,[output_dir '/' sample '_' 'Cell_Labels' '.tif']); %Saves connection map

disp('End of single cell timecourse extraction')

%% 3. Timecourse Detrending

%normalizing values; cell by cell
% detrend and normalize the calcium for all cells

nums = [1:length(CellTC(1,:))];

for cell = 1:length(CellTC(1,:))   

    current = CellTC(:,cell);
    norm_current = detrend(current,"omitnan");
    trend = current - norm_current;
    avg_detrended = mean(norm_current);
 
figure, plot(current); %plots intensity timecourse of entire islet so user can identify first responder and wave origin ranges
hold on
plot(norm_current)
plot(trend)
yline(avg_detrended)
legend("OG signal","detrended data","trend","mean of detrended", "location","northwest")
title('cell ',num2str(nums(cell)))


% save this detrended timetrace
for i = 1:length(norm_current);
    detrended_ca(i,cell) = norm_current(i);
    smooth_detrended(i,cell) = smoothdata(norm_current(i));
end 
end 


pause(3)
close all

%% Duty Cycle

% Choose starting and ending low frames; 
% determine what amplitude is considered "on"

figure, plot(detrended_ca)
figure, plot(smooth_detrended)
pause(3)

%thresh = input('Input Max Baseline Value\n'); %(x coordinate value)
%manually input what the threshold should be
%on = thresh; % correct for if avg_amp is slightly biased

for cell = 1:length(detrended_ca(1,:))
     current = detrended_ca(:,cell);
     current2 = smoothdata(current);
      thresh = mean(current);
      on = thresh; 
    % see what this looks like
     figure, 
     plot (current)
     hold on
     plot (current2)
     yline(thresh,'--','baseline/thresh')
     title('cell ',num2str(nums(cell)))
    
    % find what percent of the time the cell is "on"
    logical = find(current2 > on);
    duty_cycle(cell) = length(logical)/length(current);

 end 

 duty_cycle = duty_cycle';

%% Amplitude (as a factor of normalized trend)
cell_peaks = zeros(20,length(detrended_ca(1,:)));

 for cell = 1:length(detrended_ca(1,:))
    Mean = smooth_detrended(:,cell);
    [B,indx] = maxk(Mean,20);
    highestfreq = Mean(indx);
    avg_amplitude(cell) = mean(highestfreq);
    cell_peaks(:,cell) = indx;
 end 

%% Frequency 
 for cell = 1:length(CellTC(1,:))
   % Frequency
    Mean = CellTC(:,cell);
    s = Mean-mean(Mean); %use calcium to determine phase
    s = detrend(s);
    s = smooth(s);
    t = (0:length(s)-1); %time vector
    Ts = mean(diff(t)); %% Sampling period
    Fs = 1/Ts;  % Sampling Frequency
    Fn = Fs/2;  %Nyquist Frequency
    L  = length(s); %length of signal

    fts = fft(s-mean(s));
    course =s-mean(s);
    f = Fs*(0:(ceil(L/2)))/L;
    P2 = abs(fts/L);
    P1 = P2(1:(ceil(L/2)+1));
    P1(2:end-1) = 2*P1(2:end-1);

    figure, plot (0:(Fs/L):(Fs/2-Fs/L),P1(i,1:L/2))
    title ('FFT')

    amp_fts = abs(fts);      % Spectrum Amplitude
    phs_fts = angle(fts);
    highestfreq = f(find(P1==max(P1))); %this is the prob
    value = find(f==highestfreq);
    frequencyvec(cell,1) = f(value);
    amp(cell,1) = amp_fts(value);

    figure
    %plot(t, Mean-mean(Mean), 'r')
    plot(t,s,'r')
    hold on
    plot((t),max(P1)*sin(2*pi*f(value) *(t)), 'b')
    title('Fourier Frequency')

 end 

 pause(3)
 close all

 %% Compile Results
SingleCellData.dutycycle = duty_cycle;
SingleCellData.amplitude = avg_amplitude';
SingleCellData.frequency = frequencyvec;
SingleCellData.fft = driving_fft;
SingleCellData.smoothdata = smooth_detrended;