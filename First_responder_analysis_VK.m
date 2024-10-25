
%% THIS PROGRAM CALCULATES TIME OF THE 1-ST PHASE RESPONCE OF THE CALCIUM SIGNAL OF THE ISLET AND EACH INDIVIDUAL CELL 
%% INPUT PARAMETERS: NUMBER OF CELLS, TIME INTERVAL, CALCIUM TIMELAPSE DATA 
%% CREDIT: VIRA KRAVETS AUG, 2019

%% 1. INPUT PARAMETERS
filepath = 'C:\Users\Josh\Desktop\Vira_Lab_Folder\Simulations\HubSimulation\EASD_2019\Islet1_10percgk_VK\Islet1_leadingVKcorrectSeed1\';
calciumT = importdata([filepath '\calcium.txt']);

RandomVarsT=importdata([filepath '\RandomVars.txt']);

%ZeroCoupCell = importdata([filepath '\ZeroCoupCell.txt']);

numcells=1000;                                       % enter number of cells  

figure(1)
plot(calciumT)
st=100;                                          % starting frame
ed=725;                                          % ending frame

figure(3)
plot(calciumT(st:ed,:))
[Mx,IndMx]=maxk(calciumT(st:ed,:),1,1);             % returns max values of the Ca intensity for each cell, and corresponding index (time point)
[Mn,IndMn]=mink(calciumT(st:ed,:),1,1);

% 3. MAKING THE REFERENCE SIGNAL TO COMPARE THE SIGNAL OF INDIVIDUAL CELL'S Ca WITH THIS REFERENCE

for i=st:ed                                       % itterative index for time
currenttime = calciumT(i,:);                      % picking the row of the calcium table, corresponding to current timepoint, i
MeanIslet(i-(st-1)) = mean(currenttime);          % reference signal. Index (i-(st-1)) is here to account for times when st is not 0, otherwise indexing is wrong

end
figure(4)
plot(MeanIslet)                                   % plotting reference signal (Islet-average)

% 4. NORMALIZING Ca OF ISLET_AVERAGE AND EACH INDIVIDUAL CELL TO BE BETWEEN [0:1]
% and OBTAINING CROSS-CORRELATION OF THE REFERENCE SIGNAL (MEANISLET)WITH EACH INDIVIDUAL CELL

MeanIsletN=((MeanIslet-MeanIslet(1,1))./(mean(Mx)-MeanIslet(1,1)));  % normalizing islet-average to be between [0:1]
[HHval,HHtime] = min(abs(MeanIsletN-0.5));                   % HHtime - time at which Ca elevation of the Islet-Average reaches it's half-height; HHval - not important (equals to [normalized Ca intensity - 0.5])

for k=1:numcells                                  % itterative index for cells
    currentcell = calciumT(st:ed,k);              % picking column of the calcium table corresponding to current cell, j
    currentcellN =((currentcell-Mn(k))./(Mx(k)-Mn(k)));   % normalizing each cell to be between [0:1]
    calciumN(1:ed-st+1,k)=currentcellN;           % writing each normalized cell into array 
    
%     [c(:,k) lags]=xcov(currentcellN,MeanIsletN,'coeff');    % cross-covariance  measures the similarity between currentcell and shifted (lagged) copies of MeanIslet as a function of the lag.
%     [maxCV maxCL]=max(c);                                   % maxCV - maximum value of xcov; maxCL - index (lag) of the maxCV 
    
    [cHHval,cHHtime(:,k)] = min(abs(currentcellN-0.5));     % cHHtime - time at which Ca elevation of the k-th cell reaches it's half-height; HHval - not important
end

figure(5)
plot(calciumN)                                               % plotting normalized [0:1] Ca timelapses

[FirstVal,FirstInd]=mink(cHHtime,1000,2);                    % FirstVal - HHElevT of the 1000 (or other number) 1st-responders; FirstInd - cell numbers of the first 1000 (or other number) 1st-responders
FirstIndCol=FirstInd';                                       % converting cell numbers of the 1st responding cells to column vector
Tresp=cHHtime';                                              % converting time at which Ca elevation of the k-th cell reaches it's half-height to column vector

% 5. PLOTTING SIGNAL, XCOV, AND OUTPUTTING MAX XCOV AND CORRESPONDING TIME LAG

% figure(5)
% plot(By, 'Linewidth',2)
% hold on
% plot(A - mean(A),'Linewidth',2)
% title('Islet-Average Ca: Experimental and Low-Frequency ifft fit')                                                  % Period of the main component
% 
% 6. PLOTTING gKatp and kGlyc FOR THE 1st RESPONDERS AND FOR THE ISLET-AVERAGE.
gKatp=RandomVarsT(:,1);              % pulling Katp (potassium channel conductance) from RandomVars's column 1;
gKatpMean=mean(gKatp);

gCoup=RandomVarsT(:,2);              % pulling gCoup (coupling conductance) from RandomVars's column 2;
gCoupMean=mean(gCoup);

kGlyc=RandomVarsT(:,10);              % pulling kGlyc (glycolysis rate) from RandomVars's column 10;
kGlycMean=mean(kGlyc);

%%%%%Locating first responders and their corresponding metabolic and
%%%%%electrical parameters:
for j=FirstInd(1:100)                 % here FirstInd(1:x) includes x=max number of first-responder cells to average over
    gKatpCellsF(j,:)=gKatp(j);         % locating gKatp corresponding to the first x fastest-responding cells
    gCoupCellsF(j,:)=gCoup(j);
    kGlycCellsF(j,:)=kGlyc(j);         % locating kGlyc corresponding to the first x fastest-responding cells
end
gKatpCellsNZ=nonzeros(gKatpCellsF);    % removing all the cells not belonging to fastest [1:x]
FirstgKatpMean=mean(gKatpCellsNZ);    % finding avreage gKatp for the first [1:x] cells

gCoupCellsNZ=nonzeros(gCoupCellsF);    % removing all the cells not belonging to fastest [1:x]
FirstgCoupMean=mean(gCoupCellsNZ);   % finding avreage gCoup for the first [1:x] cells

kGlycCellsNZ=nonzeros(kGlycCellsF);    % removing all the cells not belonging to fastest [1:x]
FirstgkGlycMean=mean(kGlycCellsNZ);   % finding avreage kGlyc for the first [1:x] cells

%%%%%Locating last responders and their corresponding metabolic and
%%%%%electrical parameters:
for j=FirstInd(900:1000)                 % here FirstInd(x:numcells) includes x= number of last-responder cells to average over
    gKatpCellsL(j,:)=gKatp(j);         % locating gKatp corresponding to the first x last-responding cells
    gCoupCellsL(j,:)=gCoup(j);
    kGlycCellsL(j,:)=kGlyc(j);         % locating kGlyc corresponding to the first x last-responding cells
end
gKatpCellsNZ=nonzeros(gKatpCellsL);    % removing all the cells not belonging to slowest [1:x]
LastKatpMean=mean(gKatpCellsNZ);    % finding avreage gKatp for the last [1:x] cells

gCoupCellsNZ=nonzeros(gCoupCellsL);    % removing all the cells not belonging to slowest [1:x]
LastgCoupMean=mean(gCoupCellsNZ);    % finding avreage gKatp for the last [1:x] cells

kGlycCellsNZ=nonzeros(kGlycCellsL);    % removing all the cells not belonging to slowest [1:x]
LastgkGlycMean=mean(kGlycCellsNZ);   % finding avreage kGlyc for the last [1:x] cells

figure(6)
scatter(Tresp,gKatp);
title('gKatp, Katp conductance')  
figure(7)
scatter(Tresp,kGlyc);
title('GLYCv, Glycolysis Rate')

%%%%%%%%%
FirstRespCa=calciumT(:,FirstInd(1));          %finding the ca timelapse corresponding to 1st-responder
LastRespCa=calciumT(:,FirstInd(end));         %finding the ca timelapse corresponding to last-responder
figure(8)
plot(FirstRespCa);
title('Calcium for 1st, mean, and last -responders')
hold on
plot(mean(calciumT,2));
hold on
plot(LastRespCa)