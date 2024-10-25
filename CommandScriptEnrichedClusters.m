%% %%%% Beta Cell Cluster %%%% %%
clearvars
clear Run
 close all
clc

%% Set Path and Folder Names %%
selectBackgroundCell = 1;
useGFPasMaskBoolean = 1;

pcpath = 'F:\Clusters2\';
pcpath = 'C:\Users\dwuletj\Documents\JennsFiles\MicroscopeFiles\Clusters2\';
macpath = '/Users/jdwulet/Desktop/Clusters2/';
resultsfolder = 'results_GFPonly';

folder = {'2017_04_17', '2017_04_20'};

%% Didnt analyze the following because they look dead?
%'2017_02_09/Cluster_2G_I1_timecourse'
%Cluster_2G_I3/11G - islets looks like its dying but still has activity
%2/13'Cluster_11G_I2','Cluster_2G_I2'

names = {'Enriched_Cluster_2G_I2','Enriched_Cluster_2G_I3',...
    'Enriched_Cluster_2G_I4','Enriched_Cluster_2G_I5','Enriched_Cluster_2G_I7',...
    'Enriched_Cluster_11G_I2','Enriched_Cluster_11G_I13',...
    'Enriched_Cluster_11G_I5','Enriched_Cluster_11G_I7',...
    ...
    'Enriched_Cluster_2G_I8', 'Enriched_Cluster_2G_I12_1', ...
    'Enriched_Cluster_2G_I12_2', 'Enriched_Cluster_2G_I67_1',...
    'Enriched_Cluster_2G_I67_2','Enriched_Cluster_2G_I345_1',...
    'Enriched_Cluster_2G_I345_2','Enriched_Cluster_2G_I345_3',...
    'Enriched_Cluster_11G_I8','Enriched_Cluster_11G_I12_1',...
    'Enriched_Cluster_11G_I12_2','Enriched_Cluster_11G_I26_1', ...
    'Enriched_Cluster_11G_I26_2','Enriched_Cluster_11G_I57_1',...
    'Enriched_Cluster_11G_I57_2','Enriched_Cluster_11G_I345_1',...
    'Enriched_Cluster_11G_I345_2','Enriched_Cluster_11G_I345_3'};
foldernum = [ones(1,9), 2*ones(1, 18)];
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
%redo j=18
for j=[15,25,26,27] %
    clear filename
    whichfolder = foldernum(j);
    fullpathname = [Path folder{whichfolder} slash names{j}];
    %Try to load masks already drawn%
    try
        load([fullpathname '.mat']);
        disp('Successfully loaded .mat file');
    catch
        disp(fullpathname);
        disp('Unable to load .mat file');
    end
    
    filename.Location = [fullpathname ending];
    filename.GFPLocation =[fullpathname ending];
    
%     try
        [filename Dat] = Run(filename, selectBackgroundCell, useGFPasMaskBoolean, 1);
        
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
