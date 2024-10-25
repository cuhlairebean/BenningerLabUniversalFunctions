% %% %%%% Beta Cell Cluster %%%% %%
clearvars
clear Run_Nosilentcell
close all
clc

addpath('C:\Users\dwuletj\Documents\GitHub\UniversalCode\')
%dbstop at 432 in Run
dbstop at 445 in Run_Nosilentcell

%% Set Path and Folder Names %%
useGFPasMaskBoolean = 0; %if you want to just the gfp positive 
tryloadingvalues = 1;
samefile = 1; %GFP and calcium are in the same file
%GFPsecond = 1; %Is GFP the first or second image?

pcpath = 'F:\ScopeBackup\Confocal_Backup\BetaCellClusters\';
%macpath = '/Volumes/JaeAnnSeagateHD/ScopeBackup/Confocal_Backup/';

%%Day 21
%folder ={'2017_08_12\Calcium'}; GFPsecond = 0; %didnt analyze this day
%folder ={'2017_08_10\Calcium'}; GFPsecond = 0;
% folder = {'2017_08_07\Calcium'};
%{'2017_02_13'};
%{'2017_02_09'};

%%Day30
%folder ={'2017_09_07'}; GFPsecond = 0;
%folder ={'2017_08_30\day30\calcium'}; GFPsecond = 0;
%{'2017_08_24\day30\Calcium'};...
%{'2017_08_17\day30calcium'};

%%no silent cell
%folder ={'2018_03_26'}; GFPsecond = 0; %all good
%folder ={'2018_03_13'}; GFPsecond = 0; %all good
%folder = {'2018_03_05'}; GFPsecond = 0; %all good
%folder = {'2018_03_01'}; GFPsecond = 0; %all good
%folder = {'2018_02_26'}; GFPsecond = 0; %really good islets
%folder = {'2018_02_26\Clone2'}; GFPsecond = 0; %really good islets
%folder = {'2017_11_13'}; GFPsecond = 0; %all good %really good islets
%folder ={'2017_09_07'}; GFPsecond = 0;
%folder ={'2017_08_30\day30\calcium'}; GFPsecond = 0;
%%day21
%%%folder ={'2017_08_17\Calcium'}; GFPsecond = 0; %%dont analyze
%folder ={'2017_08_12\Calcium'}; GFPsecond = 0;
folder ={'2017_08_10\Calcium'}; GFPsecond = 0;
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
jstart = 1; %file
kstart = 1; %islet

%% Run Analysis %%
errors = 0;
tic

for i = istart:length(folder)
    disp(['i = ' num2str(i)]);
    istart = 1;
    F = dir([Path folder{i} slash '*' ending]);
    
    for j=jstart:size(F,1)
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
        analyzefile2 = contains(name,'zstack','IgnoreCase',true);%doesnt analyze  files names snap
      
        if ~analyzefile || ~analyzefile2
            R = dir([Path folder{i} slash resultsfolder slash name '*_Data.mat']);
            if x<1
                x = input('how many islets');
                disp(fullpathname);
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
                save([savedfile '_thresh2_4Data.mat'], 'filename');
                %save([savedfile '_Dat.mat'], 'Dat');
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
save([Path 'workspace' datestr(datetime('today')) '.mat'])
toc
