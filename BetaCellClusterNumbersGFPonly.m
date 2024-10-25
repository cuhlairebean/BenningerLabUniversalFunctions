%% %%%% Beta Cell Cluster %%%% %%
close all
clear all
clc
% activity2mM = [];
% activity11mM = [];
% activitypurified = [];
% 
% avgcorr2mM = [];
% avgcorr11mM = [];
% %avgcorrpurified = [];
% 
% areacorr2mM = [];
% areacorr11mM = [];
% %areacorrpurified = [];

path = '/Users/jdwulet/Desktop/Clusters/';
path = 'C:\Users\dwuletj\Documents\JennsFiles\MicroscopeFiles\Clusters\results_GFPonly\';
path = '/Volumes/JaeAnnSeagateHD/MicroscopeFiles/Clusters/';

dates = {'2017_02_09', '2017_02_09', '2017_02_09','2017_02_09','2017_02_09'...
    '2017_02_13','2017_02_13','2017_02_13','2017_02_13','2017_02_13','2017_02_13',...
    '2017_02_20_purified/','2017_02_20_purified/','2017_02_20_purified/',...
    '2017_02_20_purified/','2017_02_20_purified/','2017_02_20_purified/',...
    };
folder = '/results_GFPonly/';

results = {'Cluster_2G_I2_19Jul2017','Cluster_2G_I3_19Jul2017','Cluster_11G_I1_19Jul2017',...
    'Cluster_11G_I3_19Jul2017', 'Cluster_11G_I4_19Jul2017',...
     'Cluster_2G_I1_19Jul2017','Cluster_2G_I1_19Jul2017','Cluster_2G_I1_19Jul2017'...
    'Cluster_11G_I1_19Jul2017','Cluster_11G_I3_19Jul2017', 'Cluster_11G_I4_19Jul2017'...
     'PurifiedCluster_2G_I1_20Jul2017','PurifiedCluster_2G_I2_20Jul2017','PurifiedCluster_2G_I3_20Jul2017'...
    'PurifiedCluster_11G_I1_20Jul2017','PurifiedCluster_11G_I2_20Jul2017','PurifiedCluster_11G_I3_20Jul2017'...
    };

for i = 1:length(results)
    load([path dates{i} folder results{i} '_Data.mat'])
    load([path dates{i} folder results{i} '_Dat.mat'])
    %[path purified2{i} '.mat']
%     areacorr(j) = max([filename.RatioCorr_1;filename.RatioCorr_2;filename.RatioCorr_3]);
%     avgcorr(j) = max([filename.AvgCorr(1);filename.AvgCorr(2);filename.AvgCorr(3)]);
%     mattcorr(j) = Dat.RMCA;
%     mattactive(j) = Dat.RAA;
    
    Data.path{i,1} = path;
    Data.filename{i,1} = results{i};
%     Data.glucose{i, 1} = '11mM';
%     Data.type{i,1} = 'enriched';
    Data.date{i,1} = dates{i};
    Data.activity(i,1) = filename.TestRatioActive;
    Data.areacorr(i,1) = max([filename.RatioCorr_1;filename.RatioCorr_2;filename.RatioCorr_3]);
    Data.avgcorr(i,1) = max([filename.AvgCorr(1);filename.AvgCorr(2);filename.AvgCorr(3)]);
    Data.mattcorr(i,1) = Dat.RMCA;
    Data.mattactive(i,1) = Dat.RAA;
end

figure;
boxplot([[Data.activity([1,2,6,7,8],1);nan],Data.activity([3,4,5,9,10,11],1)],'Labels',{'2mM','11mM'})
ylim([0 1])
title('Beta Cell Clusters Activity')

figure
boxplot([[Data.areacorr([1,2,6,7,8],1);nan],Data.areacorr([3,4,5,9,10,11],1)],'Labels',{'2mM','11mM'})
ylim([0 1])
title('Beta Cell Clusters Area Correlation')

figure;
boxplot([[Data.mattactive([1,2,6,7,8],1);nan],Data.mattactive([3,4,5,9,10,11],1)],'Labels',{'2mM','11mM'})
ylim([0 1])
title('Beta Cell Clusters Activity')

figure
boxplot([[Data.mattcorr([1,2,6,7,8],1);nan],Data.mattcorr([3,4,5,9,10,11],1)],'Labels',{'2mM','11mM'})
ylim([0 1])
title('Beta Cell Clusters Area Correlation')


%% Enriched
figure;
boxplot([Data.activity([12,13,14],1),Data.activity([15,16,17],1)],'Labels',{'2mM','11mM'})
ylim([0 1])
title('Enriched Beta Cell Clusters Activity')

figure
boxplot([Data.areacorr([12,13,14],1),Data.areacorr([15,16,17],1)],'Labels',{'2mM','11mM'})
ylim([0 1])
title('Enriched Beta Cell Clusters Area Correlation')

figure;
boxplot([Data.mattactive([12,13,14],1),Data.mattactive([15,16,17],1)],'Labels',{'2mM','11mM'})
ylim([0 1])
title('Enriched Beta Cell Clusters Activity')

figure
boxplot([Data.mattcorr([12,13,14],1),Data.mattcorr([15,16,17],1)],'Labels',{'2mM','11mM'})
ylim([0 1])
title('Enriched Beta Cell Clusters Area Correlation')



% filename = [path 'GFPanalysis-' datestr(datetime('today'))];
% T = struct2table(Data);
% writetable(T,[filename '.xlsx']);

saveAllFigsToPPT([path 'boxplot-' datestr(datetime('today'))]);
