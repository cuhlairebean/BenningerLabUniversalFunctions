%% THIS PROGRAM IMPORTS CALCIUM IMAGING FILES FOR CELL-BY-CELL ANALYSIS OF CALCIUM TRACE TO IDENTIFY FIRST RESPONDER AND WAVE ORIGIN/WAVE END CELLS
%% REAL-TIME USER INPUT REQUIRED 
%% CREDIT: VIRA KRAVETS, FEB 2019; WOLFGANG SCHLEICHER, MARCH 2019; JAEANN DWULET, MARCH 2019
% Version 2

%% HOUSEKEEPING
close all
clear all
clc
%% LOADING THE CA IMAGE FILE

R = bfopen('Select Calcium Imaging File'); % Uses bfopen program to open .czi/.lsm image files

omeMeta = R{1,4};
stackSizeX = omeMeta.getPixelsSizeX(0).getValue(); % image width, pixels
stackSizeY = omeMeta.getPixelsSizeY(0).getValue(); % image height, pixels

voxelSizeXdefaultValue = omeMeta.getPixelsPhysicalSizeX(0).value();           % returns value in default unit
voxelSizeXdefaultUnit = omeMeta.getPixelsPhysicalSizeX(0).unit().getSymbol(); % returns the default unit type
voxelSizeX = omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER); % in µm
voxelSizeXdouble = voxelSizeX.doubleValue();                                  % The numeric value represented by this object after conversion to type double
voxelSizeY = omeMeta.getPixelsPhysicalSizeY(0).value(ome.units.UNITS.MICROMETER); % in µm
voxelSizeYdouble = voxelSizeY.doubleValue();                                  % The numeric value represented by this object after conversion to type double

if length(R(:,4))>1 % Pulled from calcium imaging code; think it mounts each frame into 3d array?
    pics=R(:,1);
    for i=1:length(pics)
        pics1 = pics{i};
        pics2{i} = pics1{1};
    end
    pics=pics2';
else
    pics=R{1};
    pics=pics(:,1);
end

for i=1:length(pics)
    IMG(:,:,i)=pics{i};
end
pics={};

Images=double(IMG); % converts images to double precision
RawImg=Images(:,:,1); % assigns the first frame of the video to RawImg variable

%% CODE TO MODIFY TIME COURSES
    % This block converts the 3rd dimension of array into timecourse for
    % later reference
st=1;
%ed=size(Images,3);
[xx, yy, zz] = size(Images(:,:,:));
ed=zz;
T=st:1:ed;

images=Images(:,:,st:ed);
DataOut.images=images;

%% DECLARING IMAGE PROPERTIES

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
RGB = hsv2rgb(HSV); %converts to rgb image

OGFig = figure(1);
imshow(RGB);

%%Getting rid of Background using JD code
GrayFig = figure('Name','Remove Background');
ImGray = rgb2gray(RGB); %converts image to gray
imagesc(ImGray) %displays in default colormap; easier to see cell borders
disp('Select Background Signal')
bkgrndC=createMask(imfreehand); %allows user to select background area to be removed
close(GrayFig);
BackC=ImGray.*bkgrndC; %multiplies mask by the image
BackC(~logical(BackC))=nan; %takes inverse of masked image and sets it to nan
thresh = nanmax(nanmax(BackC)); %sets threshold based on the max of values after nan's are removed
close(OGFig);

[ImgNon, ~] = RemovingAreaswNoSignal(ImGray,thresh); %Uses function from JD to remove background based on threshold
ImGray(ImgNon)=0; %sets the removed areas to 0
NoSigFig = figure('Name','Background Removed; Draw ROI around Entire Islet');
imagesc(ImGray);

%% MASKING ISLET

% User draws ROI around islet to plot whole islet timecourse
% User then determines where the start and end points for first responder
% and wave origin analyses
disp('Draw ROI around entire islet');
ROIMask_Start = imfreehand(); % User draws ROI around islet for average intensity
ROIMask_Start = createMask(ROIMask_Start);
sz_FR=ed-st+1;

StartMaskStack = images.*ROIMask_Start;
IsletTC = shiftdim(mean(mean(StartMaskStack)),2);
IsletTCfig = figure('Name','Whole Islet Time Course');
plot(IsletTC); %plots intensity timecourse of entire islet so user can identify first responder and wave origin ranges
datacursormode on % User can click on the plot to see x and y coordinates
stFR = input('Input Starting Frame for First Responder Analysis\n'); %select where first responder analysis should begin (x coordinate value)
edFR = input('Input Ending Frame for First Responder Analysis\n'); % select where first responder analysis should end (x coordinate value)
stWO = input('Input Starting Frame for Wave Origin Analysis\n'); %select where wave origin analysis should begin (x coordinate value)
edWO = input('Input Ending Frame for Wave Origin Analysis\n'); %select where wave origin analysis should end (x coordinate value)
close(IsletTCfig);
close(NoSigFig);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% User Draws ROIs around each cell within the islet
% ROIs are saved in "Mask" array and called back to throughout analyses

NoSigFig = figure('Name','Draw ROIs Around Cells of Interest');
imagesc(ImGray);

n = 1;
k = 1;
while k > 0 
    disp('Draw ROIs Around Cells of Interest')
    ROIMask = imfreehand();
    ROIMask = createMask(ROIMask);
    Mask(:,:,n) = ROIMask;
    UImp = input('Select Additonal ROI? 1=Yes, 0=No \n'); % No preallocated number of ROIs, User input required to continue ROI drawing or to stop
    if UImp == 1
        k = k+1;
        n = n+1;
    elseif UImp == 0
        k = 0;
    end
end
DataOut.Mask = Mask;
close(NoSigFig);

%% CALCULATING/PLOTTING INTENSITY OVER TIME FOR EACH ROI 
[q,r,s_mask] = size(Mask); % gets size of MASK array, specifically 's' value as it pertains to number of ROIs
sz_FR=edFR-stFR+1; % assigns a range of frames to reflect first responder analysis range
sz_WO=edWO-stWO+1; % assigns another range of frames to reflect wave origin analysis range
sz_2=ed-st; % includes all frames for exporting of entire trace

PlotHandles = zeros(1,s_mask); 
PlotLabels = cell(1,s_mask);
CellTC = zeros(sz_2,s_mask);
TC = zeros(sz_2,1);

for i = 1:s_mask
%     maskstack = repmat(Mask(:,:,i),[1 1 sz_2]);
    MaskedIMGstack = images.*Mask(:,:,i);
    for ii = 1:sz_2
        TCnoZero = MaskedIMGstack(:,:,ii);
        TCnoZero = TCnoZero(TCnoZero>0);
        TC(ii) = mean(TCnoZero);
    end
    CellTC(:,i) = TC;
%     f = polyfit(t,CellTC(:,i),2);
%     CellTC(:,i) = CellTC(:,i)./polyval(f,t);
    PlotLabels{i} = ['Cell' num2str(i)];
end
% Plotting traces for entire time course
TCFig = figure('Name','Average Intensity Over Time');
plot(CellTC);
legend(PlotLabels);
DataOut.CellTC = CellTC;

%% TRANSLATING INPUT PARAMETERS FOR VK ANALYSES
calciumT = DataOut.CellTC;
numcells = s_mask;

%% VK FIRST RESPONDER CELL ANALYSIS

[Mx,IndMx]=maxk(calciumT(stFR:edFR,:),1,1);             % returns max values of the Ca intensity for each cell, and corresponding index (time point)
[Mn,IndMn]=mink(calciumT(stFR:edFR,:),1,1);

% 3. MAKING THE REFERENCE SIGNAL TO COMPARE THE SIGNAL OF INDIVIDUAL CELL'S Ca WITH THIS REFERENCE

for i=stFR:edFR                                       % itterative index for time
    currenttime_FR = calciumT(i,:);                      % picking the row of the calcium table, corresponding to current timepoint, i
    MeanIslet_FR(i-(stFR-1)) = mean(currenttime_FR);          % reference signal. Index (i-(st-1)) is here to account for times when st is not 0, otherwise indexing is wrong
end      
% figure(4)
% plot(MeanIslet)                                   % plotting reference signal (Islet-average)

% 4. NORMALIZING Ca OF ISLET_AVERAGE AND EACH INDIVIDUAL CELL TO BE BETWEEN [0:1]
% and OBTAINING CROSS-CORRELATION OF THE REFERENCE SIGNAL (MEANISLET)WITH EACH INDIVIDUAL CELL

MeanIsletNFR=((MeanIslet_FR-MeanIslet_FR(1,1))./(mean(Mx)-MeanIslet_FR(1,1)));  % normalizing islet-average to be between [0:1]
[HHval,HHtime] = min(abs(MeanIsletNFR-0.5));                   % HHtime - time at which Ca elevation of the Islet-Average reaches it's half-height; HHval - not important (equals to [normalized Ca intensity - 0.5])

for k=1:numcells                                  % itterative index for cells
    currentcell_FR = calciumT(stFR:edFR,k);              % picking column of the calcium table corresponding to current cell, j
    currentcellN_FR =((currentcell_FR-Mn(k))./(Mx(k)-Mn(k)));   % normalizing each cell to be between [0:1]
    calciumN(1:edFR-stFR+1,k)=currentcellN_FR;           % writing each normalized cell into array 
    
%     [c(:,k) lags]=xcov(currentcellN,MeanIsletN,'coeff');    % cross-covariance  measures the similarity between currentcell and shifted (lagged) copies of MeanIslet as a function of the lag.
%     [maxCV maxCL]=max(c);                                   % maxCV - maximum value of xcov; maxCL - index (lag) of the maxCV 
    
    [cHHval,cHHtime(:,k)] = min(abs(currentcellN_FR-0.5));     % cHHtime - time at which Ca elevation of the k-th cell reaches it's half-height; HHval - not important
end

figure(5)
plot(calciumN)                                               % plotting normalized [0:1] Ca timelapses

[FirstVal,FirstInd]=mink(cHHtime,1000,2);                    % FirstVal - HHElevT of the 1000 (or other number) 1st-responders; FirstInd - cell numbers of the first 1000 (or other number) 1st-responders
FirstIndCol=FirstInd';                                       % converting cell numbers of the 1st responding cells to column vector
Tresp=cHHtime';                                              % converting time at which Ca elevation of the k-th cell reaches it's half-height to column vector

NormResp = zeros(1,numcells);
for kk = 1:numcells
    NormResp(kk) = (Tresp(kk)-min(Tresp))/(max(Tresp)-min(Tresp));
end

FirstRespCa=calciumT(:,FirstInd(1));          %finding the ca timelapse corresponding to 1st-responder
LastRespCa=calciumT(:,FirstInd(end));         %finding the ca timelapse corresponding to last-responder
figure(8)
plot(FirstRespCa);
title('Calcium; First/Last Resp. and Islet Average')
hold on
plot(mean(calciumT,2));
hold on
plot(LastRespCa)
legend('First Responder','Islet Average','Last Responder');

%% VK WAVE ORIGIN/END CELL ANALYSIS 
calciumT = DataOut.CellTC;
SR = 1;
WO_T = stWO:1:edWO;
WO_R = edWO-stWO;

FFTsizeVar = zeros(14,1);  %preallocates array for "sqaring table"
for i=5:14
    FFTsizeVar(i) = 2^(i); %calculates 2 to the power of the cycle#
    if FFTsizeVar(i)>WO_R %if the above number is greater than or equal to the number of frames
        FFTsize = FFTsizeVar(i);%The FFT size will be assigned to that calculated value
        break
    end    
end

% FFTsize=1024;                                     % By default, the FFT size is the first equal or superior power of 2 of the window size: if you have 400 fr -> use 2^9 = 512 
FHigh=20;                                         % This is for high frequency removal (low-pass filter): ex, we care oly about 6 lowest frequencies -> Ff=6 
FLow=1;                                           % This is for low frequency removal (because for Indm=2 or less, FFTsize+2-Indm will not work)

% 2. SETTING TIME AND FREQUENCY STEPS
timeVec=WO_T;                                        % 'timeVec' = start:one-time-step:end
Fs=mean(diff(timeVec));                           % mean of differences between adjacent elements of timeVec (basically, mean time steps)
Fs=1/Fs;                                          % sampling frequency, Hz
AvBin=edWO-stWO+1;                                    % provides bin number corresponding to the delta T (xcov lag) = 0. Used to normalize lags.

% 3. MAKING THE REFERENCE SIGNAL TO COMPARE THE SIGNAL OF INDIVIDUAL CELL'S CROSSCORRELATION WITH THIS REFERENCE

for i=stWO:edWO                                       % itterative index for time
    currenttime_WO = calciumT(i,:);                      % picking the row of the calcium table, corresponding to current timepoint, i
    MeanIslet_WO(i-(stWO-1)) = mean(currenttime_WO);          % reference signal. Index (i-(st-1)) is here to account for times when st is not 0, otherwise indexing is wrong
end

% 4. OBTAINING CROSS-CORRELATION OF THE REFERENCE SIGNAL (MEANISLET)WITH EACH INDIVIDUAL CELL

for j=1:numcells                                  % itterative index for cells
    currentcell_WO = calciumT(stWO:edWO,j);              % picking column of the calcium table corresponding to current cell, j
    currentcell_WO=currentcell_WO';                     % transposing column to a row
    [c(:,j),lags]=xcov(currentcell_WO,MeanIslet_WO,'coeff');      % cross-covariance  measures the similarity between currentcell and shifted (lagged) copies of MeanIslet as a function of the lag.
    [maxCV,maxCL]=max(c);                                   % maxCV - maximum value of xcov; maxCL - index (lag) of the maxCV    
end

% 6. PERFORMING FAST FOURIER TRANSFORM ON ISLET MEAN CALCIUM TIMELAPSE  
 
%Tfilter=(256-(ed-st));                                        % ??? How does filtering works for time? Does Tfilter means period filter, or ...?
 
A = MeanIslet_WO;                                                 % Reference for phase lags calculation. 
Ya=fft(A - mean(A),FFTsize);                                                                
PowerSp=abs(Ya).^2;

FFTPowerSpecFig = figure('Name','FFT Power Spectrum of Islet Avg'); %%needs to export/save
plot(PowerSp(:,1:20),'-o','Linewidth',2);
title('fft Power spectrum of the Islet-Average');
Indm1 = input('please enter first bin number (>2) for the Islet-Average power spectrum peak\n');
Indm2 = input('please enter second bin number (>2) for the Islet-Average power spectrum peak\n');
Indm3 = input('please enter third bin number (>2) for the Islet-Average power spectrum peak\n');

% 7. FILTERING THE WHOLE FFT TO KEEP ONLY LOW FREQUENCY COMPONENTS, AND DELETING DC (0-th) COMPONENT

Ya(1+FHigh:FFTsize+1-FHigh)=0;                                 % This is "filtering" in the frequency domain, via 0-padding everything after FHigh lowest bins up to FFTsize-FHigh, (+1 to account tor the 0-th DC component).
Ya(1:FLow-1)=0;                                                    % Filtering out the DC component
Aya=real(ifft(Ya));                                            % Real part of the ifft from detrended mean signal
Ay=A./(mean(A))+Aya((1:edWO-stWO+1));

%Ay=A(Tfilter:ed-Tfilter)./(mean(A(Tfilter:ed-Tfilter))+Aya((Tfilter:ed-Tfilter)));     % ??? How does filtering works for time?

% 8a. Obtaining the FFT of the Islet-Average with the main frequency component of choice (#1)
Yb1=fft(A - mean(A),FFTsize);                                  % FFT of the Islet Mean for experimental fit #1
Yb1(1+Indm1:FFTsize+1-Indm1)=0;                                % setting to 0 all bins between Indm1 on both sides of the FFT spectrum    
Yb1(1:Indm1-1)=0;                                              % setting to 0 all bins including first and up to Indm1                                 
Yb1(FFTsize-Indm1+3:FFTsize)=0;                                % setting to 0 all bins including last and down to FFTsize-Indm1 (+1 for 0th freq, +2 also for adjustment)
Ab1=real(ifft(Yb1));
By1=A./(mean(A))+Ab1((1:edWO-stWO+1));

% 8b. Obtaining the FFT of the Islet-Average with the main frequency component of choice (#2)
Yb2=fft(A - mean(A),FFTsize);                                  % FFT of the Islet Mean for experimental fit #1
Yb2(1+Indm2:FFTsize+1-Indm2)=0;                                % setting to 0 all bins between Indm1 on both sides of the FFT spectrum    
Yb2(1:Indm2-1)=0;                                              % setting to 0 all bins including first and up to Indm2                                 
Yb2(FFTsize-Indm2+3:FFTsize)=0;                                % setting to 0 all bins including last and down to FFTsize-Indm2 (+1 for 0th freq, +2 also for adjustment)
Ab2=real(ifft(Yb2));
By2=A./(mean(A))+Ab2((1:edWO-stWO+1));

% 8c. Obtaining the FFT of the Islet-Average with the main frequency component of choice (#3)
Yb3=fft(A - mean(A),FFTsize);                                  % FFT of the Islet Mean for experimental fit #1
Yb3(1+Indm3:FFTsize+1-Indm3)=0;                                % setting to 0 all bins between Indm3 on both sides of the FFT spectrum    
Yb3(1:Indm3-1)=0;                                              % setting to 0 all bins including first and up to Indm3                                 
Yb3(FFTsize-Indm3+3:FFTsize)=0;                                % setting to 0 all bins including last and down to FFTsize-Indm1 (+1 for 0th freq, +2 also for adjustment)
Ab3=real(ifft(Yb3));
By3=A./(mean(A))+Ab3((1:edWO-stWO+1));

ExpIsletFitFig = figure('Name','Experimental Islet Fits');  %%needs to export/save
subplot(3,1,1)
plot(By1)
hold on
plot(A - mean(A))
subplot(3,1,2)
plot(By2)
hold on
plot(A - mean(A))
subplot(3,1,3)
plot(By3)
hold on
plot(A - mean(A))

% 9. CHOOSING THE FREQUENCY CORRESPONDING TO THE BEST FIT OF THE EXP ISLET_AVERAGE WITH THE IFFT
Indm = input('please choose which of the 3 bin numbers fits the experimental curve best & enter it\n');
Yb=fft(A - mean(A),FFTsize);                                  % FFT of the Islet Mean for experimental fit of choice
Yb(1+Indm:FFTsize+1-Indm)=0;                                 % setting to 0 all bins between Indm on both sides of the FFT spectrum    
Yb(1:Indm-1)=0;                                              % setting to 0 all bins including first and up to Indm                                 
Yb(FFTsize-Indm+3:FFTsize)=0;                                % setting to 0 all bins including last and down to FFTsize-Indm (+1 for 0th freq, +2 also for adjustment)
Ab=real(ifft(Yb));
By=A./(mean(A))+Ab((1:edWO-stWO+1));

%plotting MeanIslet (A) on the same graph with the ifft of choice
LowFreqExpFig = figure('Name','Islet Avg Ca: Experimental and Low-Frequency');  %%needs to export/save
plot(By, 'Linewidth',2)
hold on
plot(A - mean(A),'Linewidth',2)
title('Islet-Average Ca: Experimental and Low-Frequency ifft fit')

% 10. FINDING PHASE,FREQUENCY AND PERIOD of the ISLET AVERAGE SPECTRUM (specific component) 

PhaseMax =unwrap(angle(Yb(Indm)));                             % phase of the wave component corresponding to the max of the FFT power spectrum of the Islet average
FrMax=Indm*SR./FFTsize;                                          % frequency of the main component (corrsponding to the Max of the FFT)
TMax=1./FrMax;                                                  % Period of the main component

% 11. PERFORMING FAST FOURIER TRANSFORM ON INDIVIDUAL CELLS CALCIUM TIMELAPSES

for y=1:numcells                                  % itterative index for cells
    currentcell = calciumT(stWO:edWO,y);              % picking column of the calcium table corresponding to current cell, j
    currentcell=currentcell';                     % transposing column to a row
    Yc=fft(currentcell - mean(currentcell),FFTsize); 
    
    % filtering the FFT with the LowPass filter (using the upper limit=FHigh)
    Yc(1+FHigh:FFTsize-FHigh)=0;
    Ayc=real(ifft(Yc));
    Ac=A./(mean(A))+Ayc((1:edWO-stWO+1));
    
    % determining low-frequency limit (High-pass filter) for FFT Power spectrum
    PowerSp=abs(Yc).^2;                           
    DummyFig_4 = figure('Name','Please Input High-Pass bin Number');
    subplot(2,1,1);
    plot(calciumT(stWO:edWO,y)- mean(calciumT(stWO:edWO,y)),'Linewidth',2);
    subplot(2,1,2);
    plot(PowerSp(:,1:20),'-o','Linewidth',2);
    FLow=input('Please input high-pass bin number (bins [0-FHigh] will be removed)\n');
    
    % filtering the FFT with the HighPass filter (using the lower limit=FLow)
    Yc(1:FLow-1)=0;                                                  
    Yc(FFTsize+2-FLow:FFTsize)=0;
    
    % calculating xcov of the HighPass- and LowPass- filtered FFT 
    [cSlow(:,y),lags]=xcov(Ac,Ay,'coeff');        % cross-covariance  measures the similarity between currentcell and shifted (lagged) copies of MeanIslet as a function of the lag.
    [maxCVslow,maxCLslow]=max(cSlow);
    
    [FmC,IndmC]=max(abs(Yc).^2);                  % FmC - value, IndmC - index corresponding to maximum of the FFT power spectrum for a current cell 
    Yc(1+IndmC:FFTsize-IndmC)=0; 
    
    Ayc=real(ifft(Yc));
    Ac=A./(mean(A))+Ayc((1:edWO-stWO+1));
    Acc(:,y)=Ac';
    
    PhaseMaxCshift(:,y) =PhaseMax - unwrap(angle(Yc(IndmC)));                              % phase of the wave component corresponding to the max of the FFT power spectrum of the Islet average
    FrMaxC(:,y)=IndmC*SR./FFTsize;                                          % frequency of the main component (corrsponding to the Max of the FFT)
    TMaxC(:,y)=1./FrMaxC(:,y);                                                  % Period of the main component
end

close(DummyFig_4);

ifftFig = figure('Name','ifft Corresponding to Filtered FFT');

plot(Acc(:,1:5),'Linewidth',2);
legend('cell1','cell2','cell3','cell4','cell5');
title('ifft corresponding to the filtered fft Power spectrum');

close(ifftFig);
% rewriting the xcov lags matrix to only continue correllated cells 
for y=1:numcells                                  % itterative index for cells
if  maxCVslow(1,y)<0.7                            % if the cell is not coordinated with the rest of the cells during this pulse,
        maxCLslow(1,y)=AvBin;                      % asign the lag for this cell to be 0
end
end

maxCLslowN=maxCLslow-AvBin;                       % normalized xcov lags (positive means leading cell, negative means lagging cell)

% rewriting the Phase, Frequency, and Period to only continue correllated cells 
for y=1:numcells                                  % itterative index for cells
if  maxCLslowN(1,y) == 0                             % if there is no time lag,
        PhaseMaxCshift(:,y)=0;                  % asign phase of this cell's oscillation to be equal to the phase of the islet-average
        FrMaxC(:,y)=FrMax;
        TMaxC(:,y)=TMax;
end
end

% 7. PLOTTING XCOV of FFT-filtered (Slow oscillation component) signal, AND OUTPUTTING MAX XCOV(Slow) AND CORRESPONDING TIME LAG
SlowOscFig = figure('Name','Xcov of Slow Oscillations');
plot(lags, cSlow(:,(1:5)),'Linewidth',2);
legend('cell1','cell2','cell3','cell4','cell5');
title('Coordination(xcov) of the slow oscillations');

close(SlowOscFig);

XcovLagFig = figure('Name','XCoV Lags'); %%needs to export/save first subplot
plot(maxCLslowN,'-o');
title('Xcov Lags for each cell');

NormLags = zeros(1,numcells);

for kk = 1:numcells
    NormLags(kk) = (maxCLslowN(kk) - min(maxCLslowN))/(max(maxCLslowN)-min(maxCLslowN));
end

% for kk = 1:numcells
%     maxCLslowN(kk) = maxCLslowN(kk) + abs(min(maxCLslowN));
% end
% 
% for kk = 1:numcells
%     NormLags(kk) = (maxCLslowN(kk)/max(maxCLslowN));  
% end

maxCLslowN=maxCLslow-AvBin;

%% GENERATES HEATMAPS OF RESPONDERS AND WAVE CELLS 

BWFinal = zeros(xx,yy);
for j = 1:s_mask
    DummyBW = bwperim(DataOut.Mask(:,:,j));
    BWoutline(:,:,j) = DummyBW;
    STATS = regionprops(DummyBW);
    Cent(j,:) = STATS.Centroid;
end    
CellMask = double(zeros(xx,yy));

for i = 1:numcells
    CellMask = CellMask + DataOut.Mask(:,:,i).*i;
    CellMask(find(CellMask>i)) = i;      
end

WVCellMask = CellMask;
FSCellMask = CellMask;
for i = 1:numcells
    WVCellMask(find(WVCellMask == i)) = NormLags(i)+0.0001;
    FSCellMask(find(FSCellMask == i)) = NormResp(i)+0.0001;
end

MnMxLim = round(numcells/10);

for kk = 1:MnMxLim
    [WaveMin(kk), WvMinInd(kk)] = min(NormLags);
    [WaveMax(kk), WvMaxInd(kk)] = max(NormLags);
    [RespMin(kk), RspMinInd(kk)] = min(NormResp);
    [RespMax(kk), RspMaxInd(kk)] = max(NormResp);
    NormLags(WvMinInd(kk)) = NaN;
    NormLags(WvMaxInd(kk)) = NaN;
    NormResp(RspMinInd(kk)) = NaN;
    NormResp(RspMaxInd(kk)) = NaN;
end

WaveRef = [WvMinInd, WvMaxInd];
RespRef = [RspMinInd, RspMaxInd];
WaveCells = zeros(xx,yy);
RespCells = zeros(xx,yy);

for kk = 1:length(WaveRef)
    WaveCells = or(DataOut.Mask(:,:,WaveRef(kk)),WaveCells);
    RespCells = or(DataOut.Mask(:,:,RespRef(kk)),RespCells);
    if kk == length(WaveRef)
        WavePerim = bwperim(WaveCells);
        RespPerim = bwperim(RespCells);
    end
end


Wcs = [0,0,1];
Wce = [1,0,1];
WaveR = linspace(Wcs(1),Wce(1))';
WaveG = linspace(Wcs(2),Wce(2))';
WaveB = linspace(Wcs(3),Wce(3))';
WaveMap = [WaveR,WaveG,WaveB];

Fcs = [0.8,0,0];
Fce = [0,0.8,0.1];
FstR = linspace(Fcs(1),Fce(1))';
FstG = linspace(Fcs(2),Fce(2))';
FstB = linspace(Fcs(3),Fce(3))';
FirstMap = [FstR,FstG,FstB];

pause on

ImgWave = rgb2gray(RGB);
ImgWave(WavePerim) = 255;
FinalWave = imoverlayNurin(ImgWave,WVCellMask,[0.0001,1],[0,0.8],WaveMap);
colorbar;
% WaveLimImg = imoverlay(double(FinalWave),WavePerim);
% colorbar;
ci = 1;
for c = 1:j 
    text(Cent(ci,1),Cent(ci,2),num2str(c),'Color','w');
    ci = ci+1;
end

pause on

% WaveImg = rgb2gray(RGB);
% WaveImg(WavePerim) = 255;
% WaveImg = imoverlayNurin(ImgWave,WVCellMask,[0.0001,1],[],WaveMap);
% WaveLimImg = imoverlayNurin(FinalWave,WavePerim,[],[]);
% ci = 1;
% for c = 1:j 
%     text(Cent(ci,1),Cent(ci,2),num2str(c),'Color','w');
%     ci = ci+1;
% end

ImgResp = rgb2gray(RGB);
ImgResp(RespPerim) = 255;
FinalResp = imoverlayNurin(ImgResp,FSCellMask,[0.0001,1.0001],[0,1],FirstMap);
colorbar;
ci = 1;
for c = 1:j 
    text(Cent(ci,1),Cent(ci,2),num2str(c),'Color','w');
    ci = ci+1;
end

pause on

% RespLimImg = imoverlay(FinalResp,RespPerim,'y');
% ci = 1;
% for c = 1:j 
%     text(Cent(ci,1),Cent(ci,2),num2str(c),'Color','w');
%     ci = ci+1;
% end

% Pause
%% MEASURING DISTANCE BETWEEN CELLS OF INTEREST

DistFig = figure('Name','Measure Distance Between Cells of Interest');
imshow(RGB);
ci = 1;
for c = 1:j 
    text(Cent(ci,1),Cent(ci,2),num2str(c),'Color','w');
    ci = ci+1;
end

n = 1;
kk = 1;
while kk > 0
    disp('Input IDs of Cells you wish to Link');
    CellA = num2str(input('ID of Cell A \n')); % User inputs cell number for first cell
    CellB = num2str(input('ID of Cell B \n')); % User inputs cell number for next cell
    DataOut.CellDist.Label(n,1) = cellstr([CellA 'to' CellB]); % Stores the above user inputs at labels
    disp('Draw Line linking Nuclei of Cells of Interest') % Prompts user to draw line between cells they identified previously
    DistLine = drawline(); % Allows user to draw lines linking cells
    DataOut.CellDist.Value(n,1) = pdist(DistLine.Position,'euclidean').*voxelSizeXdouble; % Calculates distance of each line in micrometers
    UImp = input('Measure Another Distance? 1=Yes, 0=No \n'); % No preallocated number of ROIs, User input required to continue ROI drawing or to stop
    if kk == UImp  
        n = n+1;
    elseif UImp == 0
        kk = 0;
    end
end

CellDist = [DataOut.CellDist.Label num2cell(DataOut.CellDist.Value)]; % Concatenates distance calculation array and label array
%% SAVING AND EXPORTING RELEVANT DATA

disp('Select Folder for Output Directory');
output_dir = uigetdir;
filename = input('Input Filename\n','s');

xlswrite([output_dir '\' filename 'Distances'],CellDist); % Exports the concatenated distance and label array as excel sheet. 
saveas(FinalWave,[output_dir '\' filename 'WaveCells' '.tif']);
saveas(FinalResp,[output_dir '\' filename 'Responders' '.tif']);
saveas(TCFig,[output_dir '\' filename 'plot' '.tif']);
saveas(XcovLagFig,[output_dir '\' filename 'XCoV_Lag_Plot' '.tif']);
saveas(FFTPowerSpecFig,[output_dir '\' filename 'FFT_Power_Spec_IslAvg' '.tif']);
saveas(ExpIsletFitFig,[output_dir '\' filename 'ExpIslet_Fits' '.tif']);
saveas(LowFreqExpFig,[output_dir '\' filename 'ExpIslet_LowFreq' '.tif']);
save([output_dir '\' filename 'masks' '.mat'],'Mask');
disp('Files Successfully Saved');