clear all
close all
clc
addpath('C:\Users\dwuletj\Documents\GitHub\UniversalCode\')
addpath('/Users/jdwulet/Documents/GitHub/UniversalCode/');


% path = 'E:\BarakBlum\';
% folder = {'KO 012\KO 012\693_UCN_3_Cre_ROBO_KO-012\';
%     'KO 015\KO 015\693_UCN_3_Cre_ROBO_KO-015\';
%     'KO 027\KO 027\613_INS_1_Cre_ROBO_KO-027\';
%     'KO 028\KO 028\613_INS_1_Cre_ROBO_KO-028\';
%     'KO 030\KO 030\613_INS_1_Cre_ROBO_KO-030\';
%     'WT 021\WT 021\611_INS_1_Cre_ROBO_WT-021\';
%     'WT 023\WT 023\611_INS_1_Cre_ROBO_WT-023\';%quite a bit of movement
%     'WT 027\WT 027\611_INS_1_Cre_ROBO_WT-027\';
%     'WT 029_2\WT 029\611_INS_1_Cre_ROBO_WT-029\';...
%     'TSeries-01042019-1223-033\TSeries-01042019-1223-033\';...
%     'TSeries-01042019-1223-034\TSeries-01042019-1223-034\';...
%     'TSeries-01042019-1223-036\TSeries-01042019-1223-036\'};


% pcpath ='G:\BarakBlum2\';
% macpath = '/Volumes/Seagate Backup Plus Drive/BarakBlum2/';
% %
% folder ={'Control 1 (ID1)/Control 1 (ID1)/TSeries-04052018-0924-035/';...
%     'Control 1 (ID1)\Control 1 (ID1)\TSeries-04052018-0924-037\';...
%     'Control 2 (ID420)\Control 2 (ID420)\TSeries-04052018-0924-037\';...
%     'Control 2 (ID420)\Control 2 (ID420)\TSeries-04052018-0924-038\';...
% 'Control 2 (ID420)\Control 2 (ID420)\TSeries-04052018-0924-042\';...
% 'Robo KO1 (ID 245)\Robo KO1 (ID 245)\TSeries-04052018-0924-030\';...
% 'Robo KO1 (ID 245)\Robo KO1 (ID 245)\TSeries-04052018-0924-031\';...
% 'Robo KO1 (ID 245)\Robo KO1 (ID 245)\TSeries-04052018-0924-032\'};

%more 1sec
% pcpath ='G:\raw_data_intravital_2021\1_sec\';
% macpath = '/Volumes/Seagate Backup Plus Drive/raw_data_intravital_2021/1_sec/';
% 
% folder ={'32', '33', '35', '36', '37', '38', '39', '40', '41', '42', '43'};
% savename = 'NewSet_1sec';

% pcpath ='G:\raw_data_intravital_2021\6_sec\';
% macpath = '/Volumes/Seagate Backup Plus Drive/raw_data_intravital_2021/6_sec/';
% 
% folder ={'20', '21', '22', '23', '24', '25', '26', '27', '29', '44'};%...
% savename = 'NewSet_6sec_final';

pcpath ='G:\raw_data_intravital_2021\10_sec\';
macpath = '/Volumes/Seagate Backup Plus Drive/raw_data_intravital_2021/10_sec/';

folder ={'1', '3', '4', '5', '6', '7', '8', '10', '11', '12', '13','14', '15', '16','18', '19','611_27'};%...
savename = 'NewSet_10sec_finalplus1';

% more 30sec
% pcpath ='G:\raw_data_intravital_2021\30_sec\';
% macpath = '/Volumes/Seagate Backup Plus Drive/raw_data_intravital_2021/30_sec/';
% 
% folder ={ '45', '46', '47', '48'};
% savename = 'NewSet_30sec_v1_final';

if exist(pcpath, 'dir') == 7
    path = pcpath;
    slash = '\';
    pc=1;
    addpath('C:\Users\dwuletj\Documents\GitHub\UniversalCode\')
    folder = regexprep(folder, '/', '\');
else
    path = macpath;
    folder = regexprep(folder, '\', '/');
    slash = '/';
    pc=0;
    addpath('/Users/jdwulet/Documents/GitHub/UniversalCode/');
end

    
m = 1;
resultsfolder = 'results';
slash = '\';

if exist(pcpath, 'dir') == 7
    path = pcpath;
    slash = '\';
    pc=1;
    addpath('C:\Users\dwuletj\Documents\GitHub\UniversalCode\')
    folder = regexprep(folder, '/', '\');
else
    path = macpath;
    folder = regexprep(folder, '\', '/');
    slash = '/';
    pc=0;
    addpath('/Users/jdwulet/Documents/GitHub/UniversalCode/');
end

for i = 1:length(folder)
    istart = 1;
    fullpathname = [path folder{i} slash resultsfolder slash];
    R = dir([fullpathname '*2021.mat']);
    
    for j=1:size(R,1) %
        jstart = 1;
        load([fullpathname R(j).name]);
        
        %for k=1:length(filename.RatioActive)
        Data.date{m, 1} = folder{i};
        Data.filename{m, 1} = R(j).name;
        %  Data.activity(m, 1) = filename.RatioActive(k);
        % Data.areacorr(m, 1) = filename.CorrPeakAmpArea(k);
        %       Data.zplane(m, 1) = k;
        try
            Data.CorrPeakAmpAreaStack(m, 1) = filename.CorrPeakAmpAreaStack;
        catch
            Data.CorrPeakAmpAreaStack(m, 1) = 0;
        end
        try
            Data.CorrPeakAmpAreaStackFilt(m, 1) = filename.CorrPeakAmpAreaStackFilt;
        catch
            Data.CorrPeakAmpAreaStackFilt(m, 1) = 0;
        end
        try
            Data.ActiveAreaStack (m, 1) = filename.ActiveAreaStack;
        catch
            Data.ActiveAreaStack (m, 1) = 0;
        end
        try
            Data.CorrActiveAreaOnly (m, 1) = (1/filename.ActiveAreaStack)*filename.CorrPeakAmpAreaStack;
        catch
            Data.CorrActiveAreaOnly (m, 1) = 0;
        end
         try
            Data.maxphase (m, 1) = filename.maxphase;
        catch
            Data.maxphase (m, 1) = 0;
        end
        
        try
            phaseMap = filename.phaseMap;
            figure;
            imagesc(phaseMap);
            colorbar
            
            maxphase = max(max(phaseMap));
            minphase = min(min(phaseMap));
            difphase = maxphase-minphase;
            Data.diffphase (m, 1) = difphase;
        catch
            Data.diffphase (m, 1) = 0;
        end
        
        
        
        m = m+1;
        % end
        clear filename
    end
end


filename = [path savename datestr(datetime('today'))];

T = struct2table(Data);
writetable(T,[filename '.xlsx']);

%saveAllFigsToPPT([path 'BarakNumbers-' datestr(datetime('today'))]);
