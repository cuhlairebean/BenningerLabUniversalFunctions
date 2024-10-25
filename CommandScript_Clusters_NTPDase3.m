% %% %%%% Beta Cell Cluster %%%% %%
clearvars
clear Run_Nosilentcell
close all
clc

addpath('C:\Users\dwuletj\Documents\GitHub\UniversalCode\')
addpath('/Users/jdwulet/Documents/GitHub/UniversalCode/')
%dbstop at 432 in Run
%dbstop at 473 in Run_Nosilentcell

%% Set Path and Folder Names %%
useGFPasMaskBoolean = 0; %if you want to just the gfp positive 
tryloadingvalues = 1;
%samefile = 0; %GFP and calcium are in the same file
%GFPsecond = 0; %Is GFP the first or second image?

% pcpath = 'F:\ScopeBackup\Confocal_Backup\BetaCellClusters\';
% macpath = '/Volumes/JaeAnnSeagateHD/alreadytransferred/';
% macpath = '/Users/jdwulet/Desktop/Clusters/';
macpath = '/Volumes/SeagateBackupPlusDrive/BetaCellClusterUCP3KO/';
pcpath = 'G:\BetaCellClusterUCP3KO\';

%folder = {'2019_08_23'};
%folder = {'2019_08_19'};
% folder = {'2019_08_12'};
% %folder = {'2019_10_04'};
% %folder = {'2019_08_28'};
% %folder = {'2019_05_23Analyzebutnotgreat'};
% folder = {'2019_02_15'};

%folder = {'2021_04_30'};
%folder = {'2021_05_04'};
folder = {'2021_05_06'};

GFPsecond = 0; 
samefile = 0; 
ending = '.czi';

resultsfolder = 'results';

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

istart = 1; %folder
jstart = 1; %file %11 is 1-8 and 17-30
kstart = 1;%islet

%% Run Analysis %%
errors = 0;
tic

for i = istart:length(folder)
    disp(['i = ' num2str(i)]);
    istart = 1;
    F = dir([Path folder{i} slash '*' ending]);
    
    for j= jstart:size(F,1) 
        disp(['j = ' num2str(j)]);
        jstart = 1;
        name = F(j).name;
        name = regexprep(name, '[.].*', '');
        disp(name);
        fullpathname = [Path folder{i} slash name];
        filename.Location = [fullpathname ending];
        filename.GFPLocation =[fullpathname ending];
        
        oLocation = pwd;
        cd([Path folder{i}])
        if (exist('results', 'dir')~=7)
            mkdir 'results';
        end
        cd(oLocation);
        x=0;
        analyzefile = contains(name,'snap','IgnoreCase',true);%doesnt analyze  files names snap
        %analyzefile2 = contains(name,'zstack','IgnoreCase',true);%doesnt analyze  files names snap
        analyzefile2 = contains(name,'GFP','IgnoreCase',true);%doesnt analyze  files names snap
        analyzefile3 = contains(name,'KCL','IgnoreCase',true);%doesnt analyze  files names snap
        if ~analyzefile && ~analyzefile2 && ~analyzefile3
            R = dir([Path folder{i} slash resultsfolder slash name '*_Dat.mat']);
           
            if tryloadingvalues
                x=length(R);
            end
            if x<1
                disp(fullpathname);
                x = input('how many islets');
                
           end
            
            for k = kstart:x
                disp(['k = ' num2str(k)]);
                kstart = 1;
                
                if tryloadingvalues
                    try
                        load([Path folder{i} slash resultsfolder slash R(k).name]);
                    catch
                        disp('Unable to load .mat file');
                    end
                end
                filename.Location = [fullpathname ending];
                filename.GFPLocation =[fullpathname ending];
                
                if useGFPasMaskBoolean
                    results = [Path folder{i} slash resultsfolder 'tdtom' slash name '_islet' num2str(k) 'tdtom'];
                    %change the name of the GFP results folders
                else
                    results = [Path folder{i} slash resultsfolder slash name '_islet' num2str(k)];
                end
                
                %[filename] = Run(filename, useGFPasMaskBoolean, samefile, GFPsecond);
                [filename] = Run_Nosilentcell(filename, useGFPasMaskBoolean, samefile, GFPsecond);
                %                 totalarea = (Dat.AAA/Dat.RAA)/filename.RatioActive;
                %                 filename.totalarea = totalarea;
                %                 filename.corrarea2_matt = Dat.AMCA/totalarea;
                %                 filename.amca_matt = Dat.AMCA;
                savedfile = [results '_' time];
                savedfile = regexprep(savedfile,'-','');
                savedfile = regexprep(savedfile,' ','_');
                %save([savedfile '_thresh2_4_Data.mat'], 'filename');
                save([savedfile '_Dat.mat'], 'filename');
                pptsavedname = [results '_' time];
                pptsavedname = regexprep(pptsavedname,'-','_');
                pptsavedname = regexprep(pptsavedname,' ','_');
                saveAllFigsToPPT(pptsavedname)
                disp(['Analysis successful for ' fullpathname]);
                
                clear filename
            end
        end
    end
end
save([Path 'UCP3KO' datestr(datetime('today')) '.mat'])
toc
