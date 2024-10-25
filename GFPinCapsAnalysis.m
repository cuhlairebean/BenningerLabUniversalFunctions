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

pcpath = 'C:\Users\dwuletj\Google Drive\LabComputer\BetaClusterResults';
% '170829_d19_gfp_4x'
% '170926_d31_4x_gfp'
macpath = '/Volumes/JaeAnnSeagateHD/ScopeBackup/Confocal_Backup/';
resultsfolder = 'results';

folder = {''};

%% Didnt analyze the following because they look dead?
ending = '.tif';

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

for i = istart:length(folder)
    istart = 1;
    F = dir([Path folder{i} slash '*' ending]);
    
    for j=jstart:size(F,1) %
        jstart = 1;
        name = F(j).name
        name = regexprep(name, '[.].*', '');
        fullpathname = [Path folder{i} slash name];
        filename.Location = [fullpathname ending];

        [filename] = GFPinCaps(filename);
        
        disp(['Analysis successful for ' fullpathname]);
        
        
        clear filename
    end
end

save([Path 'workspace' datestr(datetime('today')) '.mat'])
