%% %%%% Beta Cell Cluster %%%% %%
close all
clear all
clc

%% Set Path and Folder Names %%
pcpath = 'F:\ScopeBackup\Confocal_Backup\BetaCellClusters\';
macpath = '/Volumes/JaeAnnSeagateHD/alreadytransferred/';
resultsfolder = 'results';

%Set1 2 day 30 and 3 day 21
% folder = { '2017_08_24\day30\Calcium', '2017_08_17\day30calcium', ...
%     '2017_08_07\Calcium', '2017_02_13', '2017_02_09'};
%Set2
%day 30
folder = {'2018_03_26','2018_03_13','2018_03_05','2018_03_01','2018_02_26','2018_02_26\Clone2',...
    '2017_11_13','2017_09_07','2017_08_30\day30\calcium'};

folder = {'2019_01_19'};
  
%day 21
%folder ={'2017_08_12\Calcium','2017_08_10\Calcium'};

resultsfolder = 'results';

istart = 1;
jstart = 1;

%% Didnt analyze the following because they look dead?
ending = '.czi';

if exist(pcpath, 'dir') == 7
    Path = pcpath;
    slash = '\';
    pc=1;
else
    Path = macpath;
     folder = regexprep(folder, '\', '/');
    slash = '/';
    pc=0;
end

m = 1;
for i = istart:length(folder)
    istart = 1;
    fullpathname = [Path folder{i} slash resultsfolder slash];
    R = dir([fullpathname '*_Data.mat']);
    
    for j=jstart:size(R,1) %
        jstart = 1;
        load([fullpathname R(j).name]);
        Data.date{m, 1} = folder{i};
        Data.filename{m, 1} = R(j).name;
        Data.activity(m, 1) = filename.RatioActive;
        Data.areacorr(m, 1) = filename.CorrPeakAmpArea;
        
        %DataOut.CorrPeakAmpAreaNoGFP
        
        m = m+1;
        clear filename
    end
end
 

savefile = [Path 'BetaClusterAnalysis-' datestr(datetime('today'))];
T = struct2table(Data);
writetable(T,[savefile '.xlsx']);

save(savefile);
%saveAllFigsToPPT([macpath 'boxplot-' datestr(datetime('today'))]);