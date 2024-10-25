clear all
close all
clc
addpath('C:\Users\dwuletj\Documents\GitHub\UniversalCode\')
addpath('/Users/jdwulet/Documents/GitHub/UniversalCode/');

%% more 1sec
% pcpath ='G:\raw_data_intravital_2021\1_sec\';
% macpath = '/Volumes/Seagate Backup Plus Drive/raw_data_intravital_2021/1_sec/';
% 
% folder ={ '38', '39', '41'};
% savename = 'NewSet_1sec_phasev2_redo41';
% zstacks = 1;
% second = 1;

% % %% more 6sec
% pcpath ='G:\raw_data_intravital_2021\6_sec\';
% macpath = '/Volumes/Seagate Backup Plus Drive/raw_data_intravital_2021/6_sec/';
% 
% %folder ={'20', '23', '24', '26'};%...
% folder ={'29'};%...
% savename = 'NewSet_6sec_phasev2_num29_2';
% zstacks = 3;
% second = 6;

%% more 10sec
pcpath ='G:\raw_data_intravital_2021\10_sec\';
macpath = '/Volumes/Seagate Backup Plus Drive/raw_data_intravital_2021/10_sec/';

%folder ={'1','3', '7', '10','11', '12', '13','18'};%...
folder ={'19','611_27'};%...
savename = 'NewSet_10sec_phasev2_num61127';
zstacks = 8;
second = 10;

%% more 30sec
% pcpath ='G:\raw_data_intravital_2021\30_sec\';
% macpath = '/Volumes/Seagate Backup Plus Drive/raw_data_intravital_2021/30_sec/';

% folder ={ '45', '46', '47', '48'};
% savename = 'NewSet_30sec_phase';
% zstacks = 12;

%%analysis
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

istart = 2;
jstart = 1;

for i = istart:length(folder)
    istart = 1;
    fullpathname = [path folder{i} slash resultsfolder slash];
    R = dir([fullpathname '*2021.mat']);
%     oscillation_num = 1;
    
    j=size(R,1); %
    jstart = 1;
    name = R(j).name;
    load([fullpathname name]);
    location = filename.Location;
    location = regexprep(location, 'G:\', '/Volumes/Seagate Backup Plus Drive/');
    location = regexprep(location, '\', '/');
    
    Y=bfopen(location);
    omeMeta = Y{1, 4};
    %         stackSizeX = omeMeta.getPixelsSizeX(0).getValue(); % image width, in pixels
    %         stackSizeY = omeMeta.getPixelsSizeY(0).getValue(); % image height, in pixels
    voxelSizeX = omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER); % in µm
    voxelSizeXdouble = voxelSizeX.doubleValue(); % The numeric value represented by this object after conversion to type double
    %         voxelSizeY = omeMeta.getPixelsPhysicalSizeY(0).value(ome.units.UNITS.MICROMETER); % in µm
    %         voxelSizeYdouble = voxelSizeY.doubleValue(); % The numeric value represented by this object after conversion to type double
    %
    PixelSize = voxelSizeXdouble; % in um
    
    images=double(filename.images(:,:,:,1));
    
    IsletAvg=mean(mean(images,1),2);
    IsletAvg = IsletAvg(:);
    
    figure; plot(IsletAvg/mean(IsletAvg));
    howmanyosc = input( 'how many oscillations? ');
    
    for j=1:howmanyosc
        starttime = input( 'what starttime? ');
        endtime = input( 'what endtime? ');
        %         starttime = -1;
        %         endtime = -1;
        %try
        [dataOut] = GetPhasev2(filename.images,filename.ActMask, zstacks, PixelSize, second, starttime, endtime);
        
        for z=1:zstacks
            Data.date{m, 1} = folder{i};
            Data.filename{m, 1} = name;
            Data.PixelSize(m,1) = PixelSize;
            Data.Osc_num(m,1) = j;
            Data.zstack(m,1) = z;
            
            Data.starttime(m,1) = starttime;
            Data.endtime(m,1) = endtime;
            
            Data.timediff(m,1) = dataOut.timediff(z);
            Data.frequency(m,1) = dataOut.frequency(z) ;
            
            Data.meanspeed(m,1) = dataOut.meanspeed(z);
            Data.medspeed(m,1) = dataOut.medspeed(z);
            Data.maxspeed(m,1) = dataOut.maxspeed(z);
            Data.minspeed(m,1) = dataOut.minspeed(z);
            
            Data.meandis(m,1) = dataOut.meandis(z);
            Data.meddis(m,1) = dataOut.meddis(z);
            Data.maxdis(m,1) = dataOut.maxdis(z);
            Data.mindis(m,1) = dataOut.mindis(z);
            m=m+1;
        end
    end
    % end
    saveAllFigsToPPT([path savename name datestr(datetime('today'))]);
    clear filename
end


resultsname = [path savename datestr(datetime('today'))];

T = struct2table(Data);
writetable(T,[resultsname '.xlsx']);
