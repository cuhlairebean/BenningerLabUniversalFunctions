clear all
close all
clc

path = '/Volumes/JaeAnnSeagateHD/MicroscopeFiles/EntrainedChR2/07_05_17/Results/';

files = {'I1_R1_T0_24Jul2017_Data','I1_R1_T1_24Jul2017_Data','I1_R1_T0_24Jul2017_Data',...
    'I1_R1_T0_24Jul2017_Data','I1_R1_T0_24Jul2017_Data','I1_R1_T0_24Jul2017_Data',...
    'I1_R1_T0_24Jul2017_Data'};

% X = flipud(Dat.stimMask);
%     imagesc(X)

islet = {'I1','I2','I3'};
region = {'R1','R2','R3','R4'};
time = {'T0','T1','T2'};
count = 1;
for i = 1:length(islet)
    for r = 1:length(region)
        for t = 1:length(time)
            load([path islet{i} '_' region{r} '_' time{t} '_24Jul2017_Data.mat'])
            Data.fileid{count,t} = [islet{i} '_' region{r} '_' time{t}];
            Data.naturalFreq(count,t) = Dat.naturalFreq;
            Data.naturalPeriod(count,t) = Dat.naturalPeriod;
            Data.drivingPeriod(count,t) = Dat.drivingPeriod;
            Data.amplitudeAvg(count,t) = Dat.amplitudeAvg;
            Data.corrArea(count,t) = Dat.corrArea;
            Data.pulseArea(count,t) = Dat.pulseArea;    
        end
        count = count+1;
    end
end

filename = [path 'ChR2EntrainedAnalysis-' datestr(datetime('today'))];
 T = struct2table(Data);
 writetable(T,[filename '.xlsx']);

saveAllFigsToPPT([path 'ChR2EntrainedAnalysis-' datestr(datetime('today'))]);
