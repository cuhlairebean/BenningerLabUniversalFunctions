function [NormLags] = WOWEAnalysis(calciumT,numcells,stWO,edWO,output_dir,filename)
%WoWeAnalysis - Wave Origin and Wave End Cell Identification
%   Calculates lag times of each cell noted in calciumT input within the
%   range of frames (stWO = starting frame, edWO = ending frame). These lag
%   times are normalized and exported. 
SR = 1;
WO_T = stWO:1:edWO;
WO_R = edWO-stWO;

svst = num2str(stWO);
sved = num2str(edWO);

FFTsizeVar = zeros(14,1);  %preallocates array for "squaring table"
for i=5:14
    FFTsizeVar(i) = 2^(i); %calculates 2 to the power of the cycle#
    if FFTsizeVar(i)>=WO_R %if the above number is greater than or equal to the number of frames
        FFTsize = FFTsizeVar(i);%The FFT size will be assigned to that calculated value
        break
    end
end
if FFTsize == WO_R
    edWO = edWO-1;
    WO_R = edWO-stWO;
end
clear i
clear FFTsizeVar
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
A = MeanIslet_WO;                                                 % Reference for phase lags calculation.
Ya=fft(A - mean(A),FFTsize);
PowerSp=abs(Ya).^2;

FFTPowerSpecFig = figure('Name','FFT Power Spectrum of Islet Avg'); %%needs to export/save
plot(PowerSp(:,1:20),'-o','Linewidth',2);
title('fft Power spectrum of the Islet-Average');
Indm1 = input('Enter First Bin Number >2\n');
Indm2 = input('Enter Second Bin Number >2\n');
Indm3 = input('Enter Third Bin Number >2\n');

saveas(FFTPowerSpecFig,[output_dir '\' filename '_' 'FFT_Power_Spec_IslAvg' '_' svst '_' sved '.fig']);

close (FFTPowerSpecFig);
% 7. FILTERING THE WHOLE FFT TO KEEP ONLY LOW FREQUENCY COMPONENTS, AND DELETING DC (0-th) COMPONENT

Ya(1+FHigh:FFTsize+1-FHigh)=0;                                 % This is "filtering" in the frequency domain, via 0-padding everything after FHigh lowest bins up to FFTsize-FHigh, (+1 to account tor the 0-th DC component).
Ya(1:FLow-1)=0;                                                    % Filtering out the DC component
Aya=real(ifft(Ya));                                            % Real part of the ifft from detrended mean signal
Ay=A./(mean(A))+Aya((1:edWO-stWO+1));

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

CompFig = figure('Name','Bin Number for Best Fit');
subplot(3,1,1)
plot(By1)
xL=get(gca,'xlim');
yL=get(gca,'ylim');
text(0.99*xL(1),0.99*yL(2),num2str(Indm1),'Color','r','FontSize',14)
hold on
plot(A - mean(A))
subplot(3,1,2)
plot(By2)
xL=get(gca,'xlim');
yL=get(gca,'ylim');
text(0.99*xL(1),0.99*yL(2),num2str(Indm2),'Color','r','FontSize',14)
hold on
plot(A - mean(A))
subplot(3,1,3)
plot(By3)
xL=get(gca,'xlim');
yL=get(gca,'ylim');
text(0.99*xL(1),0.99*yL(2),num2str(Indm3),'Color','r','FontSize',14)
hold on
plot(A - mean(A))

% 9. CHOOSING THE FREQUENCY CORRESPONDING TO THE BEST FIT OF THE EXP ISLET_AVERAGE WITH THE IFFT
Indm = input ('Enter Bin Number of Best Fitting Curve\n');
saveas(CompFig,[output_dir '\' filename '_' 'Comparing_Exp_Fits' '_' svst '_' sved '.fig']);
close(CompFig);
Yb=fft(A - mean(A),FFTsize);                                 % FFT of the Islet Mean for experimental fit of choice
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

saveas(LowFreqExpFig,[output_dir '\' filename '_' 'ExpIslet_LowFreq' '_' svst '_' sved '.fig']);
close(LowFreqExpFig);

% 10. FINDING PHASE,FREQUENCY AND PERIOD of the ISLET AVERAGE SPECTRUM (specific component)

PhaseMax =unwrap(angle(Yb(Indm)));                             % phase of the wave component corresponding to the max of the FFT power spectrum of the Islet average
FrMax=Indm*SR./FFTsize;                                          % frequency of the main component (corrsponding to the Max of the FFT)
TMax=1./FrMax;                                                  % Period of the main component

% 11. PERFORMING FAST FOURIER TRANSFORM ON INDIVIDUAL CELLS CALCIUM TIMELAPSES
FLow = Indm-1;

for y=1:numcells                                  % itterative index for cells
    currentcell = calciumT(stWO:edWO,y);              % picking column of the calcium table corresponding to current cell, j
    currentcell=currentcell';                     % transposing column to a row
    Yc=fft(currentcell - mean(currentcell),FFTsize);
    
    % filtering the FFT with the LowPass filter (using the upper limit=FHigh)
    Yc(1+FHigh:FFTsize-FHigh)=0;
    Ayc=real(ifft(Yc));
    Ac=A./(mean(A))+Ayc((1:edWO-stWO+1));
    
    % determining low-frequency limit (High-pass filter) for FFT Power spectrum
%     PowerSp=abs(Yc).^2;
%     
%     PowerSp(:,21:end) = [];
%     PowerSp(:,1:2) = NaN;
% 
%     [ind,Val] = max(PowerSp);
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

% rewriting the xcov lags matrix to only continue correllated cells
%get rid of for loop using find across the board
% firstrow_maxCLslow = maxCLslow(1,:);
% firstrow_maxCLslow(find(firstrow_maxCLslow>.7))=AvBin;      
% maxCLslowtest(1,:) = firstrow_maxCLslow;

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
XcovLagFig = figure('Name','XCoV Lags'); %%needs to export/save first subplot
plot(maxCLslowN,'-o');
title('Xcov Lags for each cell');

saveas(XcovLagFig,[output_dir '\' filename '_' 'XCoV_Lag_Plot' '_' svst '_' sved '.fig']);

close(XcovLagFig);
minLag = min(maxCLslowN);
maxLag = max(maxCLslowN);
NormLags = (maxCLslowN-minLag)/(maxLag-minLag);
NormLags = NormLags';
end

