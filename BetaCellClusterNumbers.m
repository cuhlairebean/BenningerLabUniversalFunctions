%% %%%% Beta Cell Cluster %%%% %%
close all
selectBackgroundCell = 1;
activity2mM = [];
activity11mM = [];
activitypurified = [];

avgcorr2mM = [];
avgcorr11mM = [];
%avgcorrpurified = [];

areacorr2mM = [];
areacorr11mM = [];
%areacorrpurified = [];

j=1;
k=1;
macpath = '/Users/jdwulet/Desktop/Clusters/';
path = 'C:\Users\dwuletj\Documents\JennsFiles\MicroscopeFiles\Clusters\';

purified2 = {'2017_02_20_purified/PurifiedCluster_2G_I1_time';'2017_02_20_purified/PurifiedCluster_2G_I2_time';...
    '2017_02_20_purified/PurifiedCluster_2G_I3_time';'2017_02_20_purified/PurifiedCluster_2G_I4_time';};
purified11 = {'2017_02_20_purified/PurifiedCluster_11G_I1_time';'2017_02_20_purified/PurifiedCluster_11G_I2_time';...
    '2017_02_20_purified/PurifiedCluster_11G_I3_time'};

for i =1:length(purified2)
    try
    load([macpath purified2{i} '.mat'])
    [macpath purified2{i} '.mat']
    activity2mM(j) = filename.TestRatioActive;
    areacorr2mM(j) = max([filename.RatioCorr_1;filename.RatioCorr_2;filename.RatioCorr_3]);
    avgcorr2mM(j) = max([filename.AvgCorr(1);filename.AvgCorr(2);filename.AvgCorr(3)]);
    j=j+1;
    catch me
    end
end

for i =1:length(purified11)
    load([macpath purified11{i} '.mat'])
    activity11mM(k) = filename.TestRatioActive;
    areacorr11mM(k) = max([filename.RatioCorr_1;filename.RatioCorr_2;filename.RatioCorr_3]);
    avgcorr11mM(k) = max([filename.AvgCorr(1);filename.AvgCorr(2);filename.AvgCorr(3)]);
    k=k+1;
end
figure
boxplot([activity2mM',[activity11mM]'],'Labels',{'2mM','11mM'})
ylim([0 1])
title('Purified Beta Cell Clusters Activity')
figure
boxplot([areacorr2mM',[areacorr11mM]'],'Labels',{'2mM','11mM'})
ylim([0 1])
title('Purified Beta Cell Clusters Area Correlation')
figure
boxplot([avgcorr2mM',[avgcorr11mM]'],'Labels',{'2mM','11mM'})
ylim([0 1])
title('Purified Beta Cell Clusters Average Correlation')

jpur = j;
kpur = k;

name2 = {'2017_02_09/Cluster_2G_I1_timecourse';'2017_02_09/Cluster_2G_I2_timecourse';'2017_02_09/Cluster_2G_I3_timecourse';...
    '2017_02_13/Cluster_2G_I1_timecourse';'2017_02_13/Cluster_2G_I2_timecourse';'2017_02_13/Cluster_2G_I3_timecourse';'2017_02_13/Cluster_2G_I4_timecourse'};

name11 = {'2017_02_09/Cluster_11G_I1_timecourse';'2017_02_09/Cluster_11G_I3_timecourse';'2017_02_09/Cluster_11G_I4_timecourse';...
    '2017_02_13/Cluster_11G_I1_timecourse';'2017_02_13/Cluster_11G_I2_timecourse';'2017_02_13/Cluster_11G_I3_timecourse';'2017_02_13/Cluster_11G_I4_timecourse'};

for i =1:length(name2)
    load([macpath name2{i} '.mat'])
    activity2mM(j) = filename.TestRatioActive;
    areacorr2mM(j) = max([filename.RatioCorr_1;filename.RatioCorr_2;filename.RatioCorr_3]);
    avgcorr2mM(j) = max([filename.AvgCorr(1);filename.AvgCorr(2);filename.AvgCorr(3)]);
    j=j+1;
end

for i =1:length(name11)
    load([macpath name11{i} '.mat'])
    activity11mM(k) = filename.TestRatioActive;
    areacorr11mM(k) = max([filename.RatioCorr_1;filename.RatioCorr_2;filename.RatioCorr_3]);
    avgcorr11mM(k) = max([filename.AvgCorr(1);filename.AvgCorr(2);filename.AvgCorr(3)]);
    k=k+1;
end

figure;
boxplot([activity2mM',[activity11mM]'],'Labels',{'2mM','11mM'})
ylim([0 1])
title('All Beta Cell Clusters Activity')
figure
boxplot([areacorr2mM',[areacorr11mM]'],'Labels',{'2mM','11mM'})
ylim([0 1])
title('All Beta Cell Clusters Area Correlation')
figure
boxplot([avgcorr2mM',[avgcorr11mM]'],'Labels',{'2mM','11mM'})
ylim([0 1])
title('All Beta Cell Clusters Avg Correlation Across Islet')

figure;
boxplot([activity2mM(jpur:end)',activity11mM(kpur:end)'],'Labels',{'2mM','11mM'})
ylim([0 1])
title('Beta Cell Clusters Activity')
figure
boxplot([areacorr2mM(jpur:end)',areacorr11mM(kpur:end)'],'Labels',{'2mM','11mM'})
ylim([0 1])
title('Beta Cell Clusters Area Correlation')
figure
boxplot([avgcorr2mM(jpur:end)',avgcorr11mM(kpur:end)'],'Labels',{'2mM','11mM'})
ylim([0 1])
title('Beta Cell Clusters Avg Correlation Across Islet')

saveAllFigsToPPT([macpath 'boxplotppt'])

