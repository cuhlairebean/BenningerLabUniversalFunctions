% %% %%%% 3D zstack coorelation analysis %%%% %%
clearvars
close all
clc

%%add universal functions
dbstop at 185 in Run_3D_corrall

%% Set Path and Folder Names %%
tryloadingvalues = 1;
cachannel = 3;%Is Calcium the first or second image?
howmanychannel = 3;
zstacks = 12;
combinezstacks = 0;
%first channel is vasculature
%second channel is DAPI?
%third channel is calcium

%% %%%%%%%%%%%%%%%%%%%%%%%% set 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%pcpath = 'E:\BarakBlum\KO 012\KO 012\693_UCN_3_Cre_ROBO_KO-012\';
%pcpath = 'E:\BarakBlum\KO 015\KO 015\693_UCN_3_Cre_ROBO_KO-015\';
%pcpath = 'E:\BarakBlum\KO 027\KO 027\613_INS_1_Cre_ROBO_KO-027\';
%pcpath = 'E:\BarakBlum\KO 028\KO 028\613_INS_1_Cre_ROBO_KO-028\';
%pcpath = 'E:\BarakBlum\KO 030\KO 030\613_INS_1_Cre_ROBO_KO-030\';
%pcpath = 'E:\BarakBlum\WT 021\WT 021\611_INS_1_Cre_ROBO_WT-021\';
%pcpath = 'E:\BarakBlum\WT 023\WT 023\611_INS_1_Cre_ROBO_WT-023\';%quite a bit of movement
%pcpath = 'E:\BarakBlum\WT 027\WT 027\611_INS_1_Cre_ROBO_WT-027\';
%pcpath = 'E:\BarakBlum\WT 029_2\WT 029\611_INS_1_Cre_ROBO_WT-029\'; %this one isnt working
pcpath ='E:\BarakBlum\';
macpath = '/Volumes/Seagate Backup Plus Drive/BarakBlum/';
%
folder ={'KO 012\KO 012\693_UCN_3_Cre_ROBO_KO-012\';...
    'KO 015\KO 015\693_UCN_3_Cre_ROBO_KO-015\';...
    'KO 027\KO 027\613_INS_1_Cre_ROBO_KO-027\';...
    'KO 028\KO 028\613_INS_1_Cre_ROBO_KO-028\';...
    'KO 030\KO 030\613_INS_1_Cre_ROBO_KO-030\';...
    'WT 021\WT 021\611_INS_1_Cre_ROBO_WT-021\';...
    'WT 029_2\WT 029\611_INS_1_Cre_ROBO_WT-029\';...
    'WT 023\WT 023\611_INS_1_Cre_ROBO_WT-023\';...
    'WT 027\WT 027\611_INS_1_Cre_ROBO_WT-027\'};%...

folder ={'TSeries-01042019-1223-033\TSeries-01042019-1223-033\'};
% 'TSeries-01042019-1223-034\TSeries-01042019-1223-034\';...
% 'TSeries-01042019-1223-036\TSeries-01042019-1223-036\';...
% 'TSeries-01042019-1223-035\TSeries-01042019-1223-035\'}

starttime = [-1,-1,-1,-1,-1,-1,-1,30,21,-1,-1,-1];
endtime = [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1];

starttime = [-1,-1,-1,-1,-1,21];
endtime =   [-1,-1,-1,-1,-1,-1];

%'WT 023\WT 023\611_INS_1_Cre_ROBO_WT-023\';... time 30
%'WT 027\WT 027\611_INS_1_Cre_ROBO_WT-027\' time 21

ending = '.xml';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Trying to pull in matt merrins data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pcpath ='G:\For Richard Benninger\';
% macpath = '';
%
% folder ={'SLL155Ucn3CreGCaMP6s10Goscillations2point7Gto10G'};
%
% starttime = [-1];
% endtime = [-1];
% ending = '.nd2';
%
% tryloadingvalues = 1;
% cachannel = 1;%Is Calcium the first or second image?
% howmanychannel = 1;
% zstacks = 21;
% combinezstacks = 0; %merge together

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% set 2 barak %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tryloadingvalues = 1;
cachannel = 3;%Is Calcium the first or second image?
howmanychannel = 3;
zstacks = 3;
combinezstacks = 0;

pcpath ='G:\BarakBlum2\';
macpath = '/Volumes/Seagate Backup Plus Drive/BarakBlum2/';
%

%     folder ={'Control 1 (ID1)/Control 1 (ID1)/TSeries-04052018-0924-035/';...
%     'Control 1 (ID1)\Control 1 (ID1)\TSeries-04052018-0924-037\';...
%     'Control 2 (ID420)\Control 2
%     (ID420)\TSeries-04052018-0924-037\'};started at t=15
folder = {'Control 2 (ID420)\Control 2 (ID420)\TSeries-04052018-0924-038\'};
 %folder ={'Control 2 (ID420)\Control 2 (ID420)\TSeries-04052018-0924-042\'};
%     'Robo KO1 (ID 245)\Robo KO1 (ID 245)\TSeries-04052018-0924-030\';...
%     'Robo KO1 (ID 245)\Robo KO1 (ID 245)\TSeries-04052018-0924-031\';...
%     'Robo KO1 (ID 245)\Robo KO1 (ID 245)\TSeries-04052018-0924-032\';}...
    
starttime = [-1, -1, -1, -1];
endtime = [-1, -1, -1, -1];
ending = '.xml';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
resultsfolder = 'results';

if exist(pcpath, 'dir') == 7
    Path = pcpath;
    slash = '\';
    pc=1;
    addpath('C:\Users\dwuletj\Documents\GitHub\UniversalCode\')
    folder = regexprep(folder, '/', '\');
else
    Path = macpath;
    folder = regexprep(folder, '\', '/');
    slash = '/';
    pc=0;
    addpath('/Users/jdwulet/Documents/GitHub/UniversalCode/');
end
time = datestr(datetime('today'));

istart = 1; %folder
jstart = 1; %file

%% Run Analysis %%
errors = 0;
tic

for i = 1:length(folder)
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
        
        oLocation = pwd;
        cd([Path folder{i}])
        if (exist('results', 'dir')~=7)
            mkdir 'results';
        end
        cd(oLocation);
        x=0;
        
        if combinezstacks
            R = dir([Path folder{i} slash resultsfolder slash name '*combine_stacks' '.mat']);
        else
            R = dir([Path folder{i} slash resultsfolder slash name '*.mat']);
        end
        
        if tryloadingvalues
            for k = 1:size(R,1)
                try
                    analyzefile = contains(R(k).name,'stack','IgnoreCase',true);%doesnt analyze  files names snap
                catch
                    analyzefile =0 ;
                end
                if combinezstacks == analyzefile
                    try
                        load([Path folder{i} slash resultsfolder slash R(k).name]);
                    catch
                        disp('Unable to load .mat file');
                    end
                end
            end
        end
        filename.Location = [fullpathname ending];
        
        
        results = [Path folder{i} slash resultsfolder slash name];
        if combinezstacks
            [filename] = Run_3D_combinezstacks(filename,cachannel,howmanychannel,zstacks);
        else
            %[filename] = Run_3D(filename,cachannel,howmanychannel,zstacks);
            [filename] = Run_3D_corrall(filename,cachannel,howmanychannel,zstacks,starttime(i),endtime(i),ending);
            %[filename] = adjstacksfortime(filename,cachannel,howmanychannel,zstacks,starttime(i),endtime(i));
            % [filename] = GetCorrelationMapInterpolation1(filename,cachannel,howmanychannel,zstacks,starttime(i),endtime(i));
            %[filename] = getPhaseVelocityManually(filename,cachannel,howmanychannel,zstacks,starttime(i),endtime(i));
            [filename] = GetCorrelationMap(filename,cachannel,howmanychannel,zstacks,starttime(i),endtime(i));
            %             [filename] = GetPhase(filename,cachannel,howmanychannel,zstacks,starttime(i),endtime(i));
            %             difphase(:,i)=filename.difphase;
        end
        
        
        savedfile = [results '_corall' time];
        if combinezstacks
            savedfile = [savedfile 'combine_stacks'];
        end
        
        save([savedfile 'interp1.mat'], 'filename');
        saveAllFigsToPPT(savedfile)
        disp(['Analysis successful for ' fullpathname]);
        close all
        %         phase(i) = filename.maxdifphase
        %         speed(i) = filename.speed
        %         filename.difphase
        %         filename.allspeed
        %         phasemean(i) = median(filename.difphase)
        %         speedmean(i) = median(filename.allspeed)
        
        clear filename
    end
end
% %xlswrite([Path 'Phase-' datestr(datetime('today'))],difphase);
% save([Path 'corrall_workspace' datestr(datetime('today')) '.mat'])
% toc
