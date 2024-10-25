clear all
close all
clc
addpath('C:\Users\dwuletj\Documents\GitHub\UniversalCode\')

%% Change these values %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pcpath = 'E:\BarakBlum\';
macpath='';

folder = {'KO 012\KO 012\693_UCN_3_Cre_ROBO_KO-012\'};
analysisname = 'Clusters';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = 1;

resultsfolder = 'results';
if exist(pcpath, 'dir') == 7
    path = pcpath;
    slash = '\';
    pc=1;
else
    path = macpath;
    folder = regexprep(folder, '\', '/');
    slash = '/';
    pc=0;
end

for i = 1:length(folder)
    istart = 1;
    fullpathname = [path folder{i} resultsfolder slash];
    R = dir([fullpathname '*.mat']);
    
    for j=1:size(R,1) %
        jstart = 1;
        load([fullpathname R(j).name]);
        
            Data.date{m, 1} = folder{i};
            Data.filename{m, 1} = R(j).name;

            try
                Data.CorrPeakAmpAreaStack(m, 1) = filename.CorrPeakAmpAreaStack;
            catch
                Data.CorrPeakAmpAreaStack(m, 1) = 0;
            end
            try
                Data.ActiveAreaStack (m, 1) = filename.ActiveAreaStack;
            catch
                Data.ActiveAreaStack (m, 1) = 0;
            end
%             try
%                 Data.CorrActiveAreaOnly (m, 1) = (1/filename.ActiveAreaStack)*filename.CorrPeakAmpAreaStack;
%             catch
%                 Data.CorrActiveAreaOnly (m, 1) = 0;
%             end
            
            m = m+1;
       % end
        clear filename
    end
end

filename = [path analysisname datestr(datetime('today'))];
T = struct2table(Data);
writetable(T,[filename '.xlsx']);

%saveAllFigsToPPT([path 'BarakNumbers-' datestr(datetime('today'))]);
