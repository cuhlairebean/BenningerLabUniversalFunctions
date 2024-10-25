clearvars
clear Run
close all
clc

addpath('C:\Users\dwuletj\Documents\GitHub\UniversalCode\')

%% Set Path and Folder Names %%

pcpath = 'F:\ScopeBackup\Confocal_Backup\ChRod\';
macpath = '';
resultsfolder = 'results';

folder = {'2017_07_05'};

names = {'I3_R1_T0', 'I3_R1_T1', 'I3_R1_T2','I3_R2_T0', 'I3_R2_T1', 'I3_R2_T2',...
    'I3_R3_T0', 'I3_R3_T1', 'I3_R3_T2','I3_R4_T0', 'I3_R4_T1', 'I3_R4_T2'};

% {'I1_R1_T0', 'I1_R1_T1', 'I1_R1_T2','I1_R2_T0', 'I1_R2_T1', 'I1_R2_T2',...
%     'I1_R3_T0', 'I1_R3_T1', 'I1_R3_T2','I1_R4_T0', 'I1_R4_T1', 'I1_R4_T2',...
%     ...
%     'I2_R1_T0', 'I2_R1_T1', 'I2_R1_T2','I2_R2_T0', 'I2_R2_T1', 'I2_R2_T2',...
%     'I2_R3_T0', 'I2_R3_T1', 'I2_R3_T2','I2_R4_T0', 'I2_R4_T1', 'I2_R4_T2',...
%     ...
%     'I3_R1_T0', 'I3_R1_T1', 'I3_R1_T2','I3_R2_T0', 'I3_R2_T1', 'I3_R2_T2',...
%     'I3_R3_T0', 'I3_R3_T1', 'I3_R3_T2','I3_R4_T0', 'I3_R4_T1', 'I3_R4_T2'};

foldernum = 1;

ending = '.czi';


resultsfolder = 'results';

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
for j=1:length(names)
    
    clear Dat
    whichfolder = foldernum;
    fullpathname = [Path folder{whichfolder} slash names{j}];
    resultspath = [Path folder{whichfolder} slash resultsfolder slash];
    
    oLocation = pwd;
    cd([Path folder{whichfolder}])
    if (exist('results', 'dir')~=7)
        mkdir 'results';
    end
    cd(oLocation);
    
    %Try to load masks already drawn%
    %     try
    %         load([resultspath names{j} '_24Jul2017_Data' '.mat']);
    %         disp('Successfully loaded .mat file');
    %         mask = Dat.mask;
    %     catch
    %         disp(fullpathname);
    %         disp('Unable to load .mat file');
    mask=[];
    %     end
    
    Location = [fullpathname ending];
    
    %     try
    [ Dat ] = OpenFilesHighGLSM800(Location,[], mask);
    %             Data.naturalFreq(count,t) = Dat.naturalFreq;
    %             Data.naturalPeriod(count,t) = Dat.naturalPeriod;
    %             Data.drivingPeriod(count,t) = Dat.drivingPeriod;
    %             Data.amplitudeAvg(count,t) = Dat.amplitudeAvg;
    %             Data.corrArea(count,t) = Dat.corrArea;
    %             Data.pulseArea(count,t) = Dat.pulseArea;
    
    figure
    savedfile = [resultsfolder '_' time];
    savedfile = regexprep(savedfile,'-','');
    savedfile = regexprep(savedfile,' ','_');
    save([savedfile '.mat'], 'Dat');
    pptsavedname = [resultsfolder '_' time];
    pptsavedname = regexprep(pptsavedname,'-','_');
    pptsavedname = regexprep(pptsavedname,' ','_');
    saveAllFigsToPPT(pptsavedname)
    disp(['Analysis successful for ' fullpathname]);
    
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


