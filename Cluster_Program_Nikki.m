% Tom Hraha <thomas.hraha@ucdenver.edu> April 24, 2012
% Modified by JaeAnn Dwulet 
% 01/06/2017
%
% Current workable code to:
%     a. calculate area of correlation for isolated cell clusters above 
%        three different threshold values between zero and one
%     b. calculate the average phase and period from a fft
%     c. automatically save the data to an excel spreadsheet

clear all
close all
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Open File
% enter the EXACT file path and whether it is 2mM
%fileName = '/Users/jdwulet/Documents/Files from Lacie/Nikki/CalciumAnalysis/JaeAnnUpdate/Sample Data/LSM 780 Data/12wkNOD_2mM_2.lsm'; 
fileName = '/Users/jdwulet/Documents/Files from Lacie/Nikki/CalciumAnalysis/JaeAnnUpdate/Sample Data/LSM 780 Data/12wkNOD_11mM_2.lsm';
%fileName = '/Users/jdwulet/Documents/Files from Lacie/Nikki/CalciumAnalysis/JaeAnnUpdate/Sample Data/Nikon Data/NOD_3_20mM.nd2';
%fileName = '/Users/jdwulet/Documents/Files from Lacie/Nikki/CalciumAnalysis/JaeAnnUpdate/Sample Data/Nikon Data/RAG_1_20mM.nd2';
%fileName = '/Users/jdwulet/Documents/Files from Lacie/Nikki/CalciumAnalysis/JaeAnnUpdate/Sample Data/Nikon Data/RAG_2_2mM.nd2';

data.filter = 0; %set =1 for noisier data and 0 for less noisy data (Nikon)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Open File and get time slices
%clusterNumber = '01';
data.Location = fileName;
data.newmask = 1;
data = OpenFile_JD(data);
TSize = length(data.TSize);
MSize = size(data.images,1);
NSize = size(data.images,2);
% function to isolate individual cell clusters using the imfree
% function. The mask will be made from the first image in the set.
% read the first file in the folder

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% HIAnalyze
data = GetMask_JD(data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Tom's analysis
% function to: a.) calculate the cross correlation per pixel of the cluster 
% versus the average for the overall cluster to determine if they are in 
% phase, and b.) calculate the FFT of the cluster from pixel intensity 

% multiply mask stack by raw data matrix to get isolated clustera data
matrixCluster = data.images;
    
% 1 X 1 X TSize reference array of mean values
% Ref = smooth(mean(mean(matrixCluster,2),1));
% matrixClusterRef = Ref - mean(Ref);
    
%%%%choose a coordinated cluster, or choose to use whole islet mean%%%
ButtonName = questdlg('Select a cluster',...
                              'Alert','Ok','Cancel','Cancel');


celldata.Location = data.Location;
celldata.images = data.images;
celldata.filter = data.filter;
celldata.newmask = 1;
celldata = OpenFile_JD(celldata);
celldata = GetMask_JD(celldata);
cellCluster =  celldata.images;   
%cluster mean
RefCluster = smooth(mean(mean(cellCluster,2),1));
%RefCluster = RefCluster(:);
%size(RefCluster)
cellClusterRef = RefCluster - mean(RefCluster);
%whole islet mean

Ref = smooth(mean(mean(matrixCluster,2),1));
%Ref = Ref(:);
matrixClusterRef = Ref - mean(Ref);

% establish matrices for values from the for loop
matrixXCorr  = zeros(MSize, NSize, 'single');
phaseMatrix  = zeros(MSize, NSize, 'single');
periodMatrix = zeros(MSize, NSize, 'single');
%powerMatrix  = zeros(MSize, NSize, 'single');
timeMatrix   = zeros(MSize, NSize, 'single');
    
% data fit to a second order polynomial
Xz = 1:TSize;
dataFit  = polyfit(Xz,(shiftdim(squeeze(matrixClusterRef),1)),2);
dataFit2 = polyval(dataFit,Xz);

% variables for fft
L = TSize;
NFFT = 2^nextpow2(L);    % usually 1024

% cross-correlation for each pixel versus the mean of the cluster and
% an fft for phase, period and power spectrum data

for m = 1:MSize
    for n = 1:NSize
        if matrixCluster(m,n) > 0
            XCorrArray = 0;
            % matrixClusterFit = matrixCluster(m,n,:) - dataFit2;
            matrixClusterFit = 0;
            matrixClusterFit = reshape(smooth(matrixCluster(m,n,:), ...
                7),1,TSize) - dataFit2;
            [XCorrArray, lags] = xcorr((matrixClusterFit - ...
            mean(matrixClusterFit)), ...
            cellClusterRef, 'coeff');
            [coeffIndex timeIndex] = max(XCorrArray);
            timeMatrix(m,n) = timeIndex;
            matrixXCorr(m,n) = coeffIndex;
        end

        % filter values corresponding to 
        if matrixCluster(m,n) > 0 
            Yfft = 0; 
            Xindex1 = 0; Xindex2 = 0; Yindex1 = 0; Yindex2 = 0; 
            Yfft = fft((matrixClusterFit - ...
                mean(matrixClusterFit)), NFFT) / L;

            % filter out the lower frequencies
            filterVal = 3;
            Yfft(1,1:filterVal) = 0;
            Yfft(1,(1024-(filterVal -1)):1024) = 0;

            % find the mean weighted frequency and max power
            [Yindex1 Xindex1] = max(Yfft(1:(NFFT/2)));
            [Yindex2 Xindex2] = max(abs(Yfft(1:(NFFT/2))));
            powerMatrix(m,n) = Xindex2;
            %phaseMatrix(m,n) = angle(Yindex1);

            % filter for values that correspond to free space noise
            if Xindex2 > 3 && Xindex2 < 200
                periodMatrix(m,n) = (NFFT / (Xindex2 - 1));
            else 
                periodMatrix(m,n) = 0;
            end
        end
    end
end
      
% Function to calculate the area of the region that has a cross correlation
% above a threshold, generally 0.5.
threshold1  = 0.25;
threshold2  = 0.50;
threshold3  = 0.75;
threshArea1 = 0.0;
threshArea2 = 0.0;
threshArea3 = 0.0;
matrixArea  = 0.0;
clustersMat = zeros(MSize, NSize);

for m = 1:MSize
    for n = 1:NSize
        if matrixCluster(m,n) > 0
            matrixArea = matrixArea + 1;
        end
        if matrixXCorr(m,n) >= threshold1
            threshArea1 = threshArea1 + 1;
            clustersMat(m,n) = 1;
        end
        if matrixXCorr(m,n) >= threshold2
            threshArea2 = threshArea2 + 1;
            clustersMat(m,n) = 2;
        end
        if matrixXCorr(m,n) >= threshold3
            threshArea3 = threshArea3 + 1;
            clustersMat(m,n) = 3;
        end
    end
end
    
% calculate the area of Xcorr coeff above threshold
correlation1 = threshArea1 / matrixArea;
correlation2 = threshArea2 / matrixArea;
correlation3 = threshArea3 / matrixArea;
%sprintf('Cluster area in pixels is . . .')
%disp(matrixArea)
%sprintf('And the percent area of correlation above 0.25 is . . .')
corrArea1 = (correlation1*100);
%disp(corrArea1)
%sprintf('And the percent area of correlation above 0.50 is . . .')
corrArea2 = (correlation2*100);
%disp(corrArea2)
%sprintf('And the percent area of correlation above 0.75 is . . .')
corrArea3 = (correlation3*100);
%disp(corrArea3)

% mean and std of cross correlation coefficients
[ROW COL XcorrVals] = find(matrixXCorr);
meanMatrixXcorr = nanmean(nonzeros(XcorrVals));
stdMatrixXcorr = nanstd(XcorrVals);
sprintf('The mean of the cross correlation coeff is ...')
disp(meanMatrixXcorr);

% creates a 1-D array of all nonzero elements in each data matrix
[row1,col1,phaseIndex]  = find(phaseMatrix);
[row2,col2,periodIndex] = find(periodMatrix);
%[row3,col3,powerIndex]  = find(powerMatrix);

% calculate mean and std of the phase and period the ouput to command line
%meanPhase = nanmean(phaseIndex);
%stdPhase = nanstd(phaseIndex);
%sprintf('The mean and std of the phase is ...')
%disp(meanPhase); disp(stdPhase);
%meanPeriod = nanmean(nonzeros(periodIndex));
%stdPeriod = nanstd(nonzeros(periodIndex));
%sprintf('The mean and std of the period is ...')
%disp(meanPeriod); disp(stdPeriod);
%meanPower = nanmean(nonzeros(powerIndex));
%stdPower = nanstd(nonzeros(powerIndex));
%sprintf('The mean and std of the power spectrum is ...')
%disp(meanPower); disp(stdPower);

% plot the cross correlation coefficient, phase and period data per pixel
figure
imagesc(matrixXCorr)
colorbar
title('Cross Correlation Map')

% figure
% imagesc(clustersMat)
% colorbar
% title('CorrelationMap')

%figure
%imagesc(timeMatrix)
%colorbar
%title('Time Index Map')

% figure
% imagesc(phaseMatrix)
% colorbar
% title('Phase Map')

% figure
% imagesc(periodMatrix)
% colorbar
% title('Period Map') 

% Save the data to an excel spreadsheet in home directory
%d = {'Cluster Size', 'Area >0.25', 'Area >0.5', 'Area >0.75', ... 
   %  'Mean Period', 'STD Period', 'Mean Phase','STD Phase', ...
   %  'Mean Xcorr', 'STD Xcorr',...
   %  'Mean Power', 'STD Power'; ...
   %  matrixArea corrArea1 corrArea2 corrArea3 ...
   %  meanPeriod stdPeriod meanPhase stdPhase ...
   %  meanMatrixXcorr stdMatrixXcorr ...
   %  meanPower stdPower};

%xlswrite([fileName '_cluster_' clusterNumber '.xls'], d)
    