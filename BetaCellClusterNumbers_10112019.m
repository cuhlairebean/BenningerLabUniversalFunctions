%% %%%% Beta Cell Cluster %%%% %%
close all
clear all
clc

%% Set Path and Folder Names %%
pcpath = 'F:\ScopeBackup\Confocal_Backup\Clusters2\clusters\';
pcpath = 'G:\BetaCellClusterUCP3KO\';
macpath = '/Users/jdwulet/Desktop/';
resultsfolder = 'results';

% '2019_08_19'};
% '2019_08_12'};
% '2019_08_28'};

% folder = {'2019_02_15','2019_05_23Analyzebutnotgreat','2019_08_23','2019_10_04'};
% folder = {'2019_08_19', '2019_08_12','2019_08_28'};
folder = {'2021_04_30', '2021_05_04','2021_05_06'};

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
    R = dir([fullpathname '*_Dat.mat']);
    
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
 

savefile = [Path 'UCP3KO-' datestr(datetime('today'))];
T = struct2table(Data);
writetable(T,[savefile '.xlsx']);

save(savefile);
%saveAllFigsToPPT([macpath 'boxplot-' datestr(datetime('today'))]);