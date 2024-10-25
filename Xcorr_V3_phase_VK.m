
%% THIS PROGRAM CALCULATES CROSS-COVARIANCE BETWEEN THE CALCIUM SIGNAL OF THE ISLET AND EACH INDIVIDUAL CELL 
%% INPUT PARAMETERS: NUMBER OF CELLS, TIME INTERVAL, CALCIUM TIMELAPSE DATA AS TABLE, FFT SIZE, SAMPLING RATE, LOW-PASS FILTER VALUE
%% CREDIT: VIRA KRAVETS FEB, 2019

%% 1. INPUT PARAMETERS

calciumT = table2array(calcium);                  % convert 
numcells=5;                                       % enter number of cells  

st=866;                                           % starting frame
ed=1450;                                          % ending frame
SR=1;                                   % sampling rate, 1/sec
T=st:1:ed; 
FFTsize=1024;                                     % By default, the FFT size is the first equal or superior power of 2 of the window size: if you have 400 fr -> use 2^9 = 512 
FHigh=20;                                         % This is for high frequency removal (low-pass filter): ex, we care oly about 6 lowest frequencies -> Ff=6 
FLow=1;                                           % This is for low frequency removal (because for Indm=2 or less, FFTsize+2-Indm will not work)

% 2. SETTING TIME AND FREQUENCY STEPS
timeVec=T;                                        % 'timeVec' = start:one-time-step:end
Fs=mean(diff(timeVec));                           % mean of differences between adjacent elements of timeVec (basically, mean time steps)
Fs=1/Fs;                                          % sampling frequency, Hz
AvBin=ed-st+1;                                    % provides bin number corresponding to the delta T (xcov lag) = 0. Used to normalize lags.

% 3. MAKING THE REFERENCE SIGNAL TO COMPARE THE SIGNAL OF INDIVIDUAL CELL'S CROSSCORRELATION WITH THIS REFERENCE

for i=st:ed                                       % itterative index for time
currenttime = calciumT(i,:);                      % picking the row of the calcium table, corresponding to current timepoint, i
MeanIslet(i-(st-1)) = mean(currenttime);          % reference signal. Index (i-(st-1)) is here to account for times when st is not 0, otherwise indexing is wrong

end

% 4. OBTAINING CROSS-CORRELATION OF THE REFERENCE SIGNAL (MEANISLET)WITH EACH INDIVIDUAL CELL

for j=1:numcells                                  % itterative index for cells
    currentcell = calciumT(st:ed,j);              % picking column of the calcium table corresponding to current cell, j
    currentcell=currentcell';                     % transposing column to a row
   
    [c(:,j) lags]=xcov(currentcell,MeanIslet,'coeff');      % cross-covariance  measures the similarity between currentcell and shifted (lagged) copies of MeanIslet as a function of the lag.
    [maxCV maxCL]=max(c);                                   % maxCV - maximum value of xcov; maxCL - index (lag) of the maxCV 
        
end
    
% 5. PLOTTING SIGNAL, XCOV, AND OUTPUTTING MAX XCOV AND CORRESPONDING TIME LAG
    
figure(1)
% subplot(3,1,1)
plot(calciumT(st:ed,1:5),'Linewidth',2)
legend('cell1','cell2','cell3','cell4','cell5')
title('Ca timelapse, Experimental, [calciumT]')
% subplot(3,1,2)
% plot(calciumT(st:ed,6:10),'Linewidth',2)
% legend('cell6','cell7','cell8','cell9','cell10')
% subplot(3,1,3)
% plot(calciumT(st:ed,11:15),'Linewidth',2)
% legend('cell11','cell12','cell13','cell14','cell15')

% figure(2)
% subplot(3,1,1)
% plot(calciumT(st:ed,16:20),'Linewidth',2)
% legend('cell16','cell117','cell18','cell19','cell20')
% title('Ca timelapse, Experimental, [calciumT]')
% subplot(3,1,2)
% plot(calciumT(st:ed,21:25),'Linewidth',2)
% legend('cel21','cel22','cel23','cel24','cell25')
% subplot(3,1,3)
% plot(calciumT(st:ed,26:29),'Linewidth',2)
% legend('cell26','cell27','cell28','cell29','cell30')

% 6. PERFORMING FAST FOURIER TRANSFORM ON ISLET MEAN CALCIUM TIMELAPSE  
 
%Tfilter=(256-(ed-st));                                        % ??? How does filtering works for time? Does Tfilter means period filter, or ...?
 
A = MeanIslet;                                                 % Reference for phase lags calculation. 
Ya=fft(A - mean(A),FFTsize);                                                                
PowerSp=abs(Ya).^2;

figure(3)
plot(PowerSp(:,1:20),'-o','Linewidth',2)
title('fft Power spectrum of the Islet-Average')
Indm1 = input ('please enter first bin number (>2) for the Islet-Average power spectrum peak ')
Indm2 = input ('please enter second bin number (>2) for the Islet-Average power spectrum peak ')
Indm3 = input ('please enter third bin number (>2) for the Islet-Average power spectrum peak ')

% 7. FILTERING THE WHOLE FFT TO KEEP ONLY LOW FREQUENCY COMPONENTS, AND DELETING DC (0-th) COMPONENT

Ya(1+FHigh:FFTsize+1-FHigh)=0;                                 % This is "filtering" in the frequency domain, via 0-padding everything after FHigh lowest bins up to FFTsize-FHigh, (+1 to account tor the 0-th DC component).
Ya(1:FLow-1)=0;                                                % Filtering out the DC component
Aya=real(ifft(Ya));                                            % Real part of the ifft from detrended mean signal
Ay=A./(mean(A))+Aya((1:ed-st+1));

%Ay=A(Tfilter:ed-Tfilter)./(mean(A(Tfilter:ed-Tfilter))+Aya((Tfilter:ed-Tfilter)));     % ??? How does filtering works for time?

% 8a. Obtaining the FFT of the Islet-Average with the main frequency component of choice (#1)
Yb1=fft(A - mean(A),FFTsize);                                  % FFT of the Islet Mean for experimental fit #1
Yb1(1+Indm1:FFTsize+1-Indm1)=0;                                % setting to 0 all bins between Indm1 on both sides of the FFT spectrum    
Yb1(1:Indm1-1)=0;                                              % setting to 0 all bins including first and up to Indm1                                 
Yb1(FFTsize-Indm1+3:FFTsize)=0;                                % setting to 0 all bins including last and down to FFTsize-Indm1 (+1 for 0th freq, +2 also for adjustment)
Ab1=real(ifft(Yb1));
By1=A./(mean(A))+Ab1((1:ed-st+1));

% 8b. Obtaining the FFT of the Islet-Average with the main frequency component of choice (#2)
Yb2=fft(A - mean(A),FFTsize);                                  % FFT of the Islet Mean for experimental fit #1
Yb2(1+Indm2:FFTsize+1-Indm2)=0;                                % setting to 0 all bins between Indm1 on both sides of the FFT spectrum    
Yb2(1:Indm2-1)=0;                                              % setting to 0 all bins including first and up to Indm2                                 
Yb2(FFTsize-Indm2+3:FFTsize)=0;                                % setting to 0 all bins including last and down to FFTsize-Indm2 (+1 for 0th freq, +2 also for adjustment)
Ab2=real(ifft(Yb2));
By2=A./(mean(A))+Ab2((1:ed-st+1));

% 8c. Obtaining the FFT of the Islet-Average with the main frequency component of choice (#3)
Yb3=fft(A - mean(A),FFTsize);                                  % FFT of the Islet Mean for experimental fit #1
Yb3(1+Indm3:FFTsize+1-Indm3)=0;                                % setting to 0 all bins between Indm3 on both sides of the FFT spectrum    
Yb3(1:Indm3-1)=0;                                              % setting to 0 all bins including first and up to Indm3                                 
Yb3(FFTsize-Indm3+3:FFTsize)=0;                                % setting to 0 all bins including last and down to FFTsize-Indm1 (+1 for 0th freq, +2 also for adjustment)
Ab3=real(ifft(Yb3));
By3=A./(mean(A))+Ab3((1:ed-st+1));

figure(4)
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
Indm = input ('please choose which of the 3 bin numbers fits the experimental curve best & enter it')
Yb=fft(A - mean(A),FFTsize);                                  % FFT of the Islet Mean for experimental fit of choice
Yb(1+Indm:FFTsize+1-Indm)=0;                                 % setting to 0 all bins between Indm on both sides of the FFT spectrum    
Yb(1:Indm-1)=0;                                              % setting to 0 all bins including first and up to Indm                                 
Yb(FFTsize-Indm+3:FFTsize)=0;                                % setting to 0 all bins including last and down to FFTsize-Indm (+1 for 0th freq, +2 also for adjustment)
Ab=real(ifft(Yb));
By=A./(mean(A))+Ab((1:ed-st+1));

%plotting MeanIslet (A) on the same graph with the ifft of choice
figure(5)
plot(By, 'Linewidth',2)
hold on
plot(A - mean(A),'Linewidth',2)
title('Islet-Average Ca: Experimental and Low-Frequency ifft fit')

% 7. FINDING PHASE,FREQUENCY AND PERIOD of the ISLET AVERAGE SPECTRUM (specific component) 

PhaseMax =unwrap(angle(Yb(Indm)));                             % phase of the wave component corresponding to the max of the FFT power spectrum of the Islet average
FrMax=Indm*SR./FFTsize;                                          % frequency of the main component (corrsponding to the Max of the FFT)
TMax=1./FrMax;                                                  % Period of the main component

 % 6. PERFORMING FAST FOURIER TRANSFORM ON INDIVIDUAL CELLS CALCIUM TIMELAPSES

for y=1:numcells;                                  % itterative index for cells
    currentcell = calciumT(st:ed,y);              % picking column of the calcium table corresponding to current cell, j
    currentcell=currentcell';                     % transposing column to a row
    Yc=fft(currentcell - mean(currentcell),FFTsize); 
    
    % filtering the FFT with the LowPass filter (using the upper limit=FHigh)
    Yc(1+FHigh:FFTsize-FHigh)=0;
    Ayc=real(ifft(Yc));
    Ac=A./(mean(A))+Ayc((1:ed-st+1));
    
    % determining low-frequency limit (High-pass filter) for FFT Power spectrum
    PowerSp=abs(Yc).^2;                           
    figure(6)
    subplot(2,1,1)
    plot(calciumT(st:ed,y)- mean(calciumT(st:ed,y)),'Linewidth',2)
    subplot(2,1,2)
    plot(PowerSp(:,1:20),'-o','Linewidth',2)
    FLow=input('Please input high-pass bin number (bins [0-FHigh] will be removed)')
    
    % filtering the FFT with the HighPass filter (using the lower limit=FLow)
    Yc(1:FLow-1)=0;                                                  
    Yc(FFTsize+2-FLow:FFTsize)=0;
    
    % calculating xcov of the HighPass- and LowPass- filtered FFT 
    [cSlow(:,y) lags]=xcov(Ac,Ay,'coeff');        % cross-covariance  measures the similarity between currentcell and shifted (lagged) copies of MeanIslet as a function of the lag.
    [maxCVslow maxCLslow]=max(cSlow);
    
    [FmC IndmC]=max(abs(Yc).^2);                  % FmC - value, IndmC - index corresponding to maximum of the FFT power spectrum for a current cell 
    Yc(1+IndmC:FFTsize-IndmC)=0; 
    
    Ayc=real(ifft(Yc));
    Ac=A./(mean(A))+Ayc((1:ed-st+1));
    Acc(:,y)=Ac';
    
    PhaseMaxCshift(:,y) =PhaseMax - unwrap(angle(Yc(IndmC)));                              % phase of the wave component corresponding to the max of the FFT power spectrum of the Islet average
    FrMaxC(:,y)=IndmC*SR./FFTsize;                                          % frequency of the main component (corrsponding to the Max of the FFT)
    TMaxC(:,y)=1./FrMaxC(:,y);                                                  % Period of the main component
end


figure(7)
%subplot(3,1,1)
plot(Acc(:,1:5),'Linewidth',2)
legend('cell1','cell2','cell3','cell4','cell5')
title('ifft corresponding to the filtered fft Power spectrum')
% subplot(3,1,2)
% plot(Acc(:,6:10),'Linewidth',2)
% legend('cell6','cell7','cell8','cell9','cell10')
% subplot(3,1,3)
% plot(Acc(:,11:15),'Linewidth',2)
% legend('cell11','cell112','cell13','cell14','cell15')

% figure(8)
% subplot(3,1,1)
% plot(Acc(:,16:20),'Linewidth',2)
% legend('cell16','cell17','cell18','cell19','cell20')
% title('ifft corresponding to the filtered fft Power spectrum')
% subplot(3,1,2)
% plot(Acc(:,21:25),'Linewidth',2)
% legend('cel21','cel22','cel23','cel24','cell25')
% subplot(3,1,3)
% plot(Acc(:,26:29),'Linewidth',2)
% legend('cell26','cell27','cell28','cell29','cell30')

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
figure(9)
% subplot(3,1,1)
plot(lags, cSlow(:,(1:5)),'Linewidth',2)
legend('cell1','cell2','cell3','cell4','cell5')
title('Coordination(xcov) of the slow oscillations')
% subplot(3,1,2)
% plot(lags, cSlow(:,(6:10)),'Linewidth',2)
% legend('cell6','cell7','cell8','cell9','cell10')
% subplot(3,1,3)
% plot(lags, cSlow(:,(11:15)),'Linewidth',2)
% legend('cell11','cell112','cell13','cell14','cell15')

% figure(10)
% subplot(3,1,1)
% plot(lags, cSlow(:,(16:20)),'Linewidth',2)
% legend('cell16','cell17','cell18','cell19','cell20')
% title('Coordination(xcov) of the slow oscillations')
% subplot(3,1,2)
% plot(lags, cSlow(:,(21:25)),'Linewidth',2)
% legend('cel21','cell22','cell23','cell24','cell25')
% subplot(3,1,3)
% plot(lags, cSlow(:,(26:29)),'Linewidth',2)
% legend('cell26','cell27','cell28','cell29','cell30')

figure(11)
subplot(2,2,1)
plot(maxCLslowN,'-o')
title('Xcov Lags for each cell')
subplot(2,2,2)
plot(PhaseMaxCshift,'-d')
title('Phase shift compared to Islet-Average')
subplot(2,2,3)
plot(FrMaxC,'-s')
title('Frequency of the slow oscillation, Hz')
subplot(2,2,4)
plot(TMaxC,'-*')
title('Period, sec')
%%%%%%%%%%%%%%%%%%%%%%%%
%CODE BELOW NEEDS WORK%
%%%%%%%%%%%%%%%%%%%%%%%%

% T1=1;
% T2=ed - st; 
%     for y=1:numcells                                       % go through all cells points 1-by-1
%         D=calciumT(st:ed,y)';                           % 
%         Yd(y,:)=fft(D-mean(D),256);                         % fft of "detrended" D
%         Yd(Ff:256-Ff)=0;                                  % replace every frequency between Ff and 512-Ff with 0s
%         Dyd=real(ifft(Yd));                                % real part of the ifft of "filtered" signal
%         Dyd=Dyd';
%         %Dy=D(Tfilter:sz-Tfilter)./(mean(D(Tfilter:sz-Tfilter))+Dyd(Tfilter:sz-Tfilter));     % ??? How does filtering works for time?
%         Xa = xcorr(Aya(T1:T2)-mean(Aya(T1:T2)),Dyd(T1:T2)-mean(Dyd(T1:T2)));                     % xcorr between the Reference (Ay) and data (Dy)
%         %[i j k] = find(Xa==max(Xa((T2-T1)-30:(T2-T1)+30)));                                  % MAX of the xcorr. Why -30 and +30?
%         XaN  = xcorr(Aya(T1:T2)-mean(Aya(T1:T2)),Dy(T1:T2)-mean(Dy(T1:T2)),'coeff');            % ??? Normalized xcorr?
%         
%         %Xcorr phase and amp calculated directly        
%         %%Intensity(x,y) = mean(images(x,y,:),3);
%         Yd = fft(Dyd-mean(Dyd),1024);
%         Yxa = fft(Xa,1024);
%         LowF=3;                                                                              %??? what LowF means? Is it to get rid of some low frequencies?
%         [C,I] = max(abs(Yd(LowF:256)));
%         Phase2(x,y) = angle(Yxa(I+LowF-1));
%         Period(x,y)=(1024/(I+LowF-2));
%         ACComp(x,y) = sum(abs(Yd(I+LowF-2:I+LowF+2)))/sum(abs(Yd(1:256)));
%         %Check and change frequency
%         %FT xcorr to get phase of fundememntal
%     end