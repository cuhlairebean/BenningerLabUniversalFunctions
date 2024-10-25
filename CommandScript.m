%% %%%% Beta Cell Cluster %%%% %%
clearvars
clear Run
close all
clc

%% Set Path and Folder Names %%
selectBackgroundCell = 1;
useGFPasMaskBoolean = 1;

pcpath = 'F:\Clusters\';
%pcpath = 'C:\Users\dwuletj\Documents\JennsFiles\MicroscopeFiles\Clusters\';
macpath = '/Users/jdwulet/Desktop/Clusters2/';
resultsfolder = 'results_GFPonly';

folder = {'2017_02_09', '2017_02_13','2017_02_20_purified'};

%% Didnt analyze the following because they look dead?
%'2017_02_09/Cluster_2G_I1_timecourse'
%Cluster_2G_I3/11G - islets looks like its dying but still has activity
%2/13'Cluster_11G_I2','Cluster_2G_I2'

names = {'Cluster_2G_I2','Cluster_2G_I3',...
    'Cluster_11G_I1','Cluster_11G_I3','Cluster_11G_I4',...
    'Cluster_2G_I1','Cluster_2G_I3',...
    'Cluster_2G_I4','Cluster_11G_I1',...
    'Cluster_11G_I3','Cluster_11G_I4',...
    'PurifiedCluster_2G_I1','PurifiedCluster_2G_I2','PurifiedCluster_2G_I3',...
    'PurifiedCluster_11G_I1','PurifiedCluster_11G_I2','PurifiedCluster_11G_I3'};
calciumending = {'_timecourse','_timecourse','_time'};
foldernum = [1,1,1,1,1,2,2,2,2,2,2,3,3,3,3,3,3];
gfpending = {'_GFPsnap', '_GFP', '_snap_bothchannels', '_GFP', '_snap'};
gfpfolder = [1,1,1,1,1,2,2,2,3,3,3,4,4,4,5,5,5];
%%%%%%%%%%%%%%%%%
%folder = {'2017_04_20'};
%names = {'Enriched_Cluster_2G_I8', 'Enriched_Cluster_2G_I12_1', 'Enriched_Cluster_2G_I12_2','Enriched_Cluster_11G_I8','Enriched_Cluster_11G_I12_1', 'Enriched_Cluster_11G_I12_2'};%...
% 'Enriched_Cluster_2G_I67_1',  'Enriched_Cluster_2G_I67_2', 'Enriched_Cluster_2G_I345_1',  'Enriched_Cluster_2G_I345_2', ...
% 'Enriched_Cluster_2G_I345_3',  'Enriched_Cluster_11G_I8','Enriched_Cluster_11G_I12_1', 'Enriched_Cluster_11G_I12_2',...
% 'Enriched_Cluster_11G_I26_1', 'Enriched_Cluster_11G_I26_2', 'Enriched_Cluster_11G_I57_1', 'Enriched_Cluster_11G_I57_2',...
% 'Enriched_Cluster_11G_I345_1', 'Enriched_Cluster_11G_I345_2', 'Enriched_Cluster_11G_I345_3'};

ending = '.czi';

if exist(pcpath, 'dir') == 7
    Path = pcpath;
    slash = '\';
    pc=1;
else
    Path = macpath;
    slash = '/';
    pc=0;
end
time = datestr(datetime('today'));
%% Run Analysis %%
errors = 0;
tic
for j=15:length(names)
    whichfolder = foldernum(j);
    gfpfold = gfpfolder(j);
    fullpathname = [Path folder{whichfolder} slash names{j}];
    %Try to load masks already drawn%
    try
        load([fullpathname calciumending{whichfolder} '.mat']);
        disp('Successfully loaded .mat file');
    catch
        disp([fullpathname calciumending{whichfolder}]);
        disp('Unable to load .mat file');
    end
    
    filename.Location = [fullpathname calciumending{whichfolder} ending];
    filename.GFPLocation = [fullpathname gfpending{gfpfold} ending ];
    
%     try
        [filename Dat] = Run(filename, selectBackgroundCell, useGFPasMaskBoolean, 0);
        
        results = [Path folder{whichfolder} slash resultsfolder slash names{j}];
        savedfile = [results '_' time];
        savedfile = regexprep(savedfile,'-','');
        savedfile = regexprep(savedfile,' ','_');
        save([savedfile '_Data.mat'], 'filename');
        save([savedfile '_Dat.mat'], 'Dat');
        
        pptsavedname = [results '_' time];
        pptsavedname = regexprep(pptsavedname,'-','_');
        pptsavedname = regexprep(pptsavedname,' ','_');
        saveAllFigsToPPT(pptsavedname)
        
        disp(['Analysis successful for ' fullpathname]);
        
        %         ActiveAreafilename.ActiveArea
        %         filename.Core
        %         filename
        
        %         filename = [macpath 'analysis-' datestr(datetime('today'))];
        %         T = struct2table(Data);
        %         writetable(T,[filename '.xlsx']);
        
%     catch me
%         errors = errors+1;
%         errormat{errors,1} = fullpathname;
%         disp(['Analysis Error for ' fullpathname]);
%         errormat{errors,2} = me.message;
%         errormat{errors,3} = me.stack;
%         errormat{errors,4} = me.cause;
%         disp(me.message);
%         disp(me.stack);
%         disp(me.cause);
%     end
end
% save([Path 'workspace' datestr(datetime('today')) '.mat'])
toc
