%% %%%%Dissociated Beta Cell Cluster %%%% %%
clearvars
clear Run
close all
clc

%% Set Path and Folder Names %%
pcpath = 'C:\Users\dwuletj\Documents\JennsFiles\MicroscopeFiles\DissociatedClusters\';
macpath = '/Users/jdwulet/Desktop/DissociatedClusters/';

if exist(pcpath, 'dir') == 7
    Path = pcpath;
    slash = '\';
    pc=1;
else
    Path = macpath;
    slash = '/';
    pc=0;
end
folder = {'2017_01_31','2017_02_21'};

%%%%%%%%%%%%%%%%%
listname = {'HumanClusters_Area2003', 'HumanClusters_Area3003', 'HumanClusters_Area4003' };%, 'HumanClusters_Area5003', 'HumanClusters_Area5006', 'HumanClusters_Area6003', 'HumanClusters_Area6009'};
ins = 'C1-';
Ca = 'C3-';
BF = 'C2-';
ending = '.nd2';
foldernum = 1;
% 
% listname = {'Cluster9', 'Cluster10', 'Cluster11' };
% ins = 'ins_';
% Ca = 'ca_';
% BF = 'ca_';
% ending = '.czi'; %%%%%%%%%%%%%%%%%%%%%%%%check this%%%%%%%%%%%%%%%%%%%%%%%
% foldernum = 2;
%% Run Analysis %%

tic

for i =1:length(listname)
    name = listname{i};
    fullpathname = [Path folder{foldernum} slash];
    filename.Location = [fullpathname Ca name ending];
    BrightField.Location = [fullpathname BF name ending];
    GFP.Location = [fullpathname ins name ending];
    try
        load([fullpathname name '1.mat']);
    catch
    end
    filename.Location = [fullpathname Ca name ending];
    BrightField.Location = [fullpathname BF name ending];
    GFP.Location = [fullpathname ins name ending];
    
    filename = MakeMaps(BrightField, filename);
    savedfile = [fullpathname name '2.mat'];
    save(savedfile, 'filename')
    
    filename = DissasociatedCells(filename, GFP);
    savedfile = [fullpathname name '2.mat'];
    save(savedfile, 'filename')
    
    saveAllFigsToPPT([fullpathname name '2ppt'])
    clear filename
end
toc
%% %%%%%%%%%%%%%%% calculations
activity = [];
GFPintensity = [];
GFPpercent = [];
avgactivity = [];
cellnum = [];
avgpeakamp = [];
for i =1:length(listname)
    name = listname{i};
    fullpathname = [Path folder{foldernum} slash];
    
    load([fullpathname name '2.mat']);
    avgpeakamp = [avgpeakamp, filename.AveragePeakAmp];
    cellnum = [cellnum, 1:length(filename.NewRatioActive)]
    activity = [activity,filename.NewRatioActive];
    GFPintensity = [GFPintensity, filename.GFPintensityminusbg];
    
end
figure
plot(filename.cellTC(:,12))
hold on
plot(filename.cellTC(:,13))
hold on
plot(filename.cellTC(:,16))
hold on
plot(filename.cellTC(:,2))
hold on
legend({'cell1', 'cell2', 'cell3', 'cell4'}, 'Location','northeastoutside')

title('Timecourses of Selected Cells')