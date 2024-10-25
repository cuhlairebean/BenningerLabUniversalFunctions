%% THIS PROGRAM CALCULATES CROSS-COVARIANCE BETWEEN THE CALCIUM SIGNAL OF THE ISLET AND EACH INDIVIDUAL CELL
%% INPUT PARAMETERS: NUMBER OF CELLS, TIME INTERVAL, CALCIUM TIMELAPSE DATA AS TABLE

%% CREDIT: VIRA KRAVETS FEB, 2019/JAEANN DWULET
clear all
close all
clc

%addpath('C:\Users\dwuletj\Documents\GitHub\UniversalCode\');     %put your universal code here;

%25% variation in kglc
% path = 'F:\Simulations_AfterKeystone\Vira\Islet1_25percgk\';
% saveworkspacename = 'Islet1WaveInitiators_25perc';

%10% variation in kglc
path = 'F:\Simulations_AfterKeystone\Vira\Redo2\Islet1\';             %?insert the path to the folder with simulations downloaded from the Supercomputer
%path = 'C:\Users\Josh\Desktop\Vira_Lab_Folder\Simulations\HubSimulation\'; %main folder

saveworkspacename = 'Islet1WaveInitiators_test';
pptsavedname = saveworkspacename;

%load(saveworkspacename);

%%
%mainpath = [path 'WT_G11_Coup120\'];   %WT coupling                               % ...?insert subfolder name inside the Simulation folder
mainpath = [path 'Orig_120coupling\'];

calciumfile = [mainpath 'calcium.txt'];                               % calls for the calcium.txt file inside the subfolder
randomvars= [mainpath 'RandomVars.txt'];                              % name of the Randomvars.txt file inside the subfolder

calcium = importdata(calciumfile);
save('wildtype_8G_calcium2','calcium')                               % what is reffered to as "WT_11g_calcium2"?
 
figure
plot(calcium);
title('Full Time Course');

st=2000;                                          % enter start time for full time course for overall frequency of islet
ed=4000;                                          % enter end time for above
cashort1 = calcium(st:ed, :);

[freqvec_full, ~] = findavgfreq(cashort1) %finds avg freq

shortst=2000;                                          % enter start time that you will determine lagging and leading cells
shorted=4000;                                          % enter end time
cashort = calcium(shortst:shorted, :);

figure
plot(cashort);
title('Find lagging and leading on this timecourse');

%origtime = st:ed;
figure
xq = 1:.005:(ed-st+1);
vq1 = interp1(1:(ed-st+1),cashort,xq);
plot(1:(ed-st+1),cashort(:,1:100:1000),'o',xq,vq1(:,1:10),':.');
title('Linear Interpolation');

clear cashort mainpath calciumfile calcium

%% 1. CODE TO SET TCs
calciumT = vq1;                           % enter number of cells

% 2. MAKING THE REFERENCE SIGNAL TO COMPARE THE SIGNAL OF INDIVIDUAL CELL'S CROSSCORRELATION WITH THIS REFERENCE
MeanIslet= mean(calciumT,2)';       % reference signal. Index (i-(st-1)) is here to account for times when st is not 0, otherwise indexing is wrong

% 3. OBTAINING CROSS-CORRELATION OF THE REFERENCE SIGNAL (MEANISLET)WITH EACH INDIVIDUAL CELL
numcells=size(calciumT,2);
for j=1:numcells   % itterative index for cells
    
    currentcell = calciumT(:,j)';
    %refcell = calciumT(:,1)'; % picking column of the calcium table corresponding to current cell, j
    
    [c(:,j) lags]=xcov(currentcell-mean(currentcell),MeanIslet-mean(MeanIslet),'coeff');      % cross-covariance  measures the similarity between currentcell and shifted (lagged) copies of MeanIslet as a function of the lag      % cross-covariance  measures the similarity between currentcell and shifted (lagged) copies of MeanIslet as a function of the lag.
    [maxCV maxCL]=max(c);
   
end

% 4. PLOTTING SIGNAL, XCOV, AND OUTPUTTING MAX XCOV AND CORRESPONDING TIME LAG


% figure
% plot(calciumT(:,:),'Linewidth',2)
%legend('cell1','cell2','cell3','cell4','cell5','cell6','cell7','cell8','cell9','cell10')
figure
plot(lags, c,'Linewidth',2)
title('cross correlation')
%legend('cell1','cell2','cell3','cell4','cell5','cell6','cell7','cell8','cell9','cell10')
figure
plot(maxCL-mean(maxCL),'-o')
title('normalized lags')

phasevecsort = sort(maxCL);
minvalue = phasevecsort(100); %i got 101 cells
leadingcells = find(maxCL <= minvalue);
leadingcellsind = find(maxCL <= minvalue)'-1;%This is the values you pull for the c++ code
notleadingcells = find(maxCL > minvalue); 

maxvalue = phasevecsort(length(phasevecsort)-100);
laggingcells = find(maxCL >= maxvalue);
laggingcellsind = find(maxCL >= maxvalue)'-1;
notlaggingcells = find(maxCL < maxvalue); 

% clearvars -except leadingcells leadingcellsind notleadingcells laggingcells laggingcellsind notlaggingcells

%% get couping 0 frequency vs gk
clear cashort mainpath calciumfile calcium
mainpath = [path 'Orig_0coupling\']; %0 coupling

calciumfile = [mainpath 'calcium.txt'];
randomvars= [mainpath 'RandomVars.txt'];

calcium = importdata(calciumfile);
randomVars = importdata(randomvars);

save('Coup0_11G_calcium','calcium')
figure
plot(calcium);
title('No coupling');

% st=2000;                                          % enter start time full time course
% ed=4000; 

cashort = calcium(st:ed, :);

%1st column gKATP
%10 column glyc
%4 column pSERCA
gKATP = randomVars(:,1);
pSERCA = randomVars(:,4);
glyc = randomVars(:,10);

avggKATP = mean(gKATP)
avggKATP_lead = mean(gKATP(leadingcells))
avggKATP_lag = mean(gKATP(laggingcells))

avgpSERCA = mean(pSERCA)
avgpSERCA_lead = mean(pSERCA(leadingcells))
avgpSERCA_lag = mean(pSERCA(laggingcells))

avggglyc = mean(glyc)
avgglyc_lead = mean(glyc(leadingcells))
avgglyc_lag = mean(glyc(laggingcells))

[avggfreq_total, freqvec_0coupling] = findavgfreq(cashort);
avgfreq_lead = mean(freqvec_0coupling(leadingcells))
avgfreq_lag = mean(freqvec_0coupling(laggingcells));

figure
plot(glyc, freqvec_0coupling,'o')
lsline
mdl = fitlm(glyc, freqvec_0coupling');
r_glyc_squared = mdl.Rsquared.Adjusted;
r_glyc = corrcoef(glyc, freqvec_0coupling');
xlabel('glyc')
ylabel('freq')
title('glyc vs freq')

figure
plot(gKATP, freqvec_0coupling,'o')
lsline
mdl = fitlm(gKATP, freqvec_0coupling');
r_gkatp_squared = mdl.Rsquared.Adjusted;
r_gkatp = corrcoef(gKATP, freqvec_0coupling');
xlabel('gKATP')
ylabel('freq')
title('gKATP vs freq')

figure
plot(pSERCA, freqvec_0coupling,'o')
lsline
mdl = fitlm(pSERCA, freqvec_0coupling');
r_pSERCA_squared = mdl.Rsquared.Adjusted;
r_pSERCA = corrcoef(pSERCA, freqvec_0coupling');
xlabel('pSERCA')
ylabel('freq')
title('pSERCA vs freq')

% %% Remove Leading cells
% clear cashort mainpath calciumfile calcium
% mainpath = [path 'Islet1_leading\']; %islet without leading cells
% calciumfile = [mainpath 'calcium.txt'];
% 
% calcium = importdata(calciumfile);
% save('Islet1_LeadingRemoved','calcium')
% 
% figure
% plot(calcium);
% title('Islet leading Uncoupled-all cells') %looks uncoupled
% 
% figure
% plot(calcium(:, leadingcells));
% title('Only Leading Uncoupled')
% 
% figure
% plot(calcium(:, notleadingcells));
% title('Non leading Coupled')
% 
% % st=2000;                                          % enter start time full time course
% % ed=4000;                                          % enter end time
% cashort = calcium(st:end, :);
% %[freqvec_leading, ~] = findavgfreq(cashort);
% [freqvec_leadnotincludedincalc, ~] = findavgfreq(cashort(:,notleadingcells))
% 
% 
% %% Remove Lagging cells
% clear cashort mainpath calciumfile calcium
% mainpath = [path 'Islet1_lagging\'];
% calciumfile = [mainpath 'calcium.txt'];
% 
% calcium = importdata(calciumfile);
% save('Islet1_LaggingRemoved','calcium')
% 
% figure
% plot(calcium);
% title('Islet lagging uncoupled-all cells')
% 
% figure
% plot(calcium(:,laggingcells));
% title('Only Lagging Uncoupled')
% 
% figure
% plot(calcium(:,notlaggingcells));
% title('Non lagging Coupled')
% 
% % st=2000;                                          % enter start time
% % ed=4000;                                          % enter end time
% cashort = calcium(st:end, :);
% %[freqvec_lagging, ~] = findavgfreq(cashort);
% [freqvec_lagnotincludedincalc, ~] = findavgfreq(cashort(:,notlaggingcells))
% 
% %% Remove Random cells
% clear cashort mainpath calciumfile calcium
% mainpath = [path 'Islet1_random\'];
% calciumfile = [mainpath 'calcium.txt'];
% 
% %p = randperm(1000,101)-1; %c++ starts at 0 and this is the jth random
% %cells that are uncoupled
% 
% RandomCells = importdata([path 'Islet1_random\ZeroCoupCell.txt']); %this is the txt i used on the supercomputer to tell the sim which cells to uncouple
% RandomCells = (RandomCells+1)'; %because c++ starts at 0 index
% allcellsindex = zeros(1,1000);
% allcellsindex(RandomCells)=1;
% nonrandomindex = find(~allcellsindex);
% 
% calcium = importdata(calciumfile);
% save('Islet1_RandomRemoved','calcium')
% 
% figure
% plot(calcium);
% title('Random All Cells')
% 
% figure
% plot(calcium(:,RandomCells));
% title('Random Uncoupled Cells')
% 
% figure
% plot(calcium(:,nonrandomindex));
% title('Random Sim - Coupled cells')
% % 
% % st=2000;                                          % enter start time
% % ed=4000;                                          % enter end time
% 
% cashort = calcium(st:end, :);
%                                        % enter end time
% %[freqvec_random, ~] = findavgfreq(cashort);
% [freqvec_randomnotincludedincalc, ~] = findavgfreq(cashort(:,nonrandomindex))
% 
% saveAllFigsToPPT(pptsavedname)
% clear cashort mainpath calciumfile calcium c currentcell calciumT lags MeanIslet vq1 xq randomVars cashort1
% save(saveworkspacename)

