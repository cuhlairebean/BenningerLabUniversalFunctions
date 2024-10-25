% %% %%%% Beta Cell Cluster %%%% %%
clearvars
clear Run
close all
clc

%% Set Path and Folder Names %%
selectBackgroundCell = 1;
useGFPasMaskBoolean = 0;
tryloadingvalues = 0;
samefile = 1; %GFP and calcium are in the same file

pcpath = 'F:\ScopeBackup\Confocal_Backup\';

macpath = '/Volumes/JaeAnnSeagateHD/ScopeBackup/Confocal_Backup/';
resultsfolder = 'results';

folder = {'2017_11_16\'}
{'2017_08_24\day30\Calcium', '2017_08_30\day30\calcium'};
    %'2017_10_25\R0714'};

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

istart = 1;
jstart = 1;
kstart = 1;

%% Run Analysis %%
errors = 0;
tic
%redo j=18
for i = istart:1%length(folder)
    i
    istart = 1;
    F = dir([Path folder{i} slash '*' ending]);
    
    for j=jstart:1%size(F,1) %
        jstart = 1;
        name = F(j).name
        name = regexprep(name, '[.].*', '');
        name = 'R2056_11G_I5' %'R2056_11G_I5' %test 2mM
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
            kstart = 1;
            if tryloadingvalues
                load([Path folder{i} slash resultsfolder slash R(k).name])
            end
            filename.Location = [fullpathname ending];
            filename.GFPLocation =[fullpathname ending];
            results = [Path folder{i} slash resultsfolder slash name '_islet' num2str(k)];
            
            [filename Dat] = Run(filename, selectBackgroundCell, useGFPasMaskBoolean, samefile);
            totalarea = (Dat.AAA/Dat.RAA)/filename.RatioActive;
            filename.totalarea = totalarea;
            filename.corrarea2_matt = Dat.AMCA/totalarea;
            filename.amca_matt = Dat.AMCA;
             
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
            
            clear filename
        end
    end
end
save([Path 'workspace' datestr(datetime('today')) '.mat'])
toc
