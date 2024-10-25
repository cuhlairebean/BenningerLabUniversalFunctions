% %% %%%% Beta Cell Cluster %%%% %%
clearvars
clear Run
close all
clc

%% Set Path and Folder Names %%
selectBackgroundCell = 1; % pretty much always on the confocal
useGFPasMaskBoolean = 0; % uses the GFP or tdtomato as a mask
tryloadingvalues = 0;
samefile = 1; %GFP and calcium are in the same file

pcpath = 'F:\ScopeBackup\Confocal_Backup\BetaCellClusters\';%path before the folder, put a slash at the end of the path

macpath = '';
resultsfolder = 'results';

folder = {'2018_04_05'};%folders within path, seperate with commas, dont put a slash at the end

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

time = datestr(datetime('today'));
time = regexprep(time,'-','');
time = regexprep(time,' ','_');

istart = 1;%folder
jstart = 1;%file
kstart = 1;%islet

%% Run Analysis %%
errors = 0;
tic

for i = istart:length(folder)
    i
    istart = 1;
    F = dir([Path folder{i} slash '*' ending]);
    
    for j=jstart:size(F,1)
        j
        jstart = 1;
        name = F(j).name
        name = regexprep(name, '[.].*', '');
        fullpathname = [Path folder{i} slash name];
        filename.Location = [fullpathname ending];
        filename.GFPLocation =[fullpathname ending];
        
        if tryloadingvalues
            R = dir([Path folder{i} slash resultsfolder slash name '*_Data.mat']);
            x = size(R,1);
            %             disp('Successfully loaded .mat file');
        else
            x = input('how many islets');
            %             disp(fullpathname);
            %             disp('Unable to load .mat file');
        end
        
        for k = kstart:x
            k
            kstart = 1;
            if tryloadingvalues
                load([Path folder{i} slash resultsfolder slash R(k).name])
            end
            filename.Location = [fullpathname ending];
            filename.GFPLocation =[fullpathname ending];
            results = [Path folder{i} slash resultsfolder slash name '_islet' num2str(k)];
            
            [filename, Dat] = Run(filename, selectBackgroundCell, useGFPasMaskBoolean, samefile);
            totalarea = filename.Area;
            activearea = filename.ActiveArea;
            %filename.corrarea_matt = Dat.AMCA/totalarea;
            if activearea==0
                    filename.corrarea_matt = 0;
                    filename.corrarea_activeonly = 0;
            else
                filename.corrarea_matt = Dat.AMCA/totalarea;
                filename.corrarea_activeonly = Dat.AMCA/activearea;
                
            end
            filename.amca_matt = Dat.AMCA;
            
            savedfile = [results '_' time];
            %savedfile = regexprep(savedfile,'-','');
            %savedfile = regexprep(savedfile,' ','_');
            save([savedfile '_Data.mat'], 'filename');
           % save([savedfile '_Dat.mat'], 'Dat');
            pptsavedname = [results '_' time];
            %pptsavedname = regexprep(pptsavedname,'-','_');
            %pptsavedname = regexprep(pptsavedname,' ','_');
            saveAllFigsToPPT(pptsavedname)
            disp(['Analysis successful for ' fullpathname]);            
            
            clear filename
        end
    end
end
save([Path 'workspace' datestr(datetime('today')) '.mat'])
toc
