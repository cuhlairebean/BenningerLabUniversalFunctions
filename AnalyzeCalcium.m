% %% %%%% Beta Cell Cluster %%%% %%
clearvars
clear Run_Nosilentcell
clear CoordinatedAreas
close all

clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
<<<<<<< Updated upstream
%% Set Path and Folder Names %% 
addpath('C:\Users\levittcl\Documents\GitHub\UniversalCode\')

Path = 'E:\2021_10_12\'; %mainpath to folder - add slash
% folder = {'2021_09_16_Ctrl'}; %folder with images in it
folder = {'T1','T0'}; %folder with images in it - could be a list of folders
savename = 'Chr2_10_12'; %Name of saved files
redobackground = 1; %yes = 1 or no=0

istart = 1; %folder to start on
jstart = 3; %file to start on
=======
%% Set Path and Folder Names %%
addpath('/Users/levittcl/Documents/UniversalCode')

Path = '/Volumes/Claire Hard Drive/DATA/gCAMP - Calcium Imaging/2021_12_07/I1-controlcont_11mMKCl.czi'; %mainpath to folder, always add slash at the end
folder = {'Matlab Analysis'}; %folder with images in it
savename = 'islet2 11mM'; %Name of saved files
redobackground = 1; %yes = 1 or no=0

istart = 1; %folder to start on
jstart = 1; %file to start on
>>>>>>> Stashed changes
kstart = 1; %islet to start on

slash = '\'; %slash - change if using a mac- '/' or pc- '\'
samefile = 0; % is calcium *with* another channel? 1-yes 0-no
cafirst = 0; %is calcium the first channel? 1-yes 0-no, if only 1 channel doesnt matter what it equals
ending = '.czi'; %what is the type of microscope file

<<<<<<< Updated upstream
threshhold = 2; %threshhold for activity. 
% Keep threshhold the same for entire project. good thresholds are between 1.75-2.5 
=======
threshhold = 1.8; %threshhold for activity. %1.9 OR 1.8
% Keep threshhold the same for entire project. good thresholds are between 1.75-2. 
>>>>>>> Stashed changes
% Threshholds might change for different dyes, microscope, etc.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

resultsfolder = 'results';
time = datestr(datetime('today'));

tryloadingvalues = 0;

%% Run Analysis %%
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
        analyzefile2 = contains(name,'YFP','IgnoreCase',true);%doesnt analyze  files names YFP
        analyzefile3 = contains(name,'KCL','IgnoreCase',true);%doesnt analyze  files names KCL
        if ~analyzefile && ~analyzefile2 && ~analyzefile3
            R = dir([Path folder{i} slash resultsfolder slash name '*_Dat.mat']);
           
            if tryloadingvalues
                x=length(R);
            end
%             if x<1
                disp(fullpathname);
                x = input('how many islets');   
%            end
            
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
                
                [filename] = Run_Nosilentcell(filename, samefile, cafirst, redobackground, threshhold);
                
                results = [Path folder{i} slash resultsfolder slash name '_islet' num2str(k)];
                savedfile = [results '_' time];
                save([savedfile '_Dat.mat'], 'filename');
                pptsavedname = [results '_' time];
%                 pptsavedname = [name '_islet' num2str(k)]; %change this back
                saveAllFigsToPPT(pptsavedname)
                disp(['Analysis successful for ' fullpathname]);
                
                clear filename
            end
        end
    end
end
% save([Path savename datestr(datetime('today')) '.mat'])

%% %%%% Compile Numbers %%%% %%

istart = 1;
jstart = 1;

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
        %Data.plateaufraction(m, 1) = filename.plateaufraction;
        %calculating the mean and distribution of the oscillation period and plateau fraction, 
        m = m+1;
        clear filename
    end
end
 
<<<<<<< Updated upstream
savefile = [Path savename datestr(datetime('today'))];
T = struct2table(Data);
writetable(T,[savefile '.xlsx']);
save(savefile);
=======
% savefile = [Path slash savename datestr(datetime('today'))];
% T = struct2table(Data);
% writetable(T,[savefile '.xlsx']);
% save(savefile);
>>>>>>> Stashed changes

