clear all
close all
clc

%% 1 second
% catc = readtable('/Users/jdwulet/Google Drive/barakphase/TimeTrace1sec38.csv');
% catc = readtable('/Users/jdwulet/Google Drive/barakphase/TimeTrace1sec39.csv');
% catc = readtable('/Users/jdwulet/Google Drive/barakphase/TimeTrace1sec41.csv');
    %% 38
    %timecourse = timecourse(1:ceil(length(timecourse)/2));
    %% 39
    %timecourse = timecourse(200:end);
    %% 41
%     timecourse = timecourse(69:end);
% starttime = [-1, 200, 69];
% endtime = [ 'half the total', -1,-1];
% second = 1;
% whichz = 1;
% z=1;

%% 6 second
% catc = readtable('/Users/jdwulet/Google Drive/barakphase/TimeTrace6sec20z1.csv');
% second = 6;
% whichz = 1;
% z=3;
% starttime = -1;
% endtime = -1;

% catc = readtable('/Users/jdwulet/Google Drive/barakphase/TimeTrace6sec23z1.csv');
% second = 6;
% whichz = 1;
% z=3;
% starttime = -1;
% endtime = -1;

% catc = readtable('/Users/jdwulet/Google Drive/barakphase/TimeTrace6sec23z3.csv');
% second = 6;
% whichz = 3;
% z=3;
% starttime = -1;
% endtime = -1;

% catc = readtable('/Users/jdwulet/Google Drive/barakphase/TimeTrace6sec24z2.csv');
% second = 6;
% whichz = 2;
% z=3;
% starttime = -1;
% endtime = 55;

% catc = readtable('/Users/jdwulet/Google Drive/barakphase/TimeTrace6sec26z3.csv');
% second = 6;
% whichz = 3;
% z=3;
% starttime = 20;
% endtime = -1;

% catc = readtable('/Users/jdwulet/Google Drive/barakphase/TimeTrace10sec1z6.csv');
% second = 10;
% whichz = 6;
% z=8;
% starttime = -1;
% endtime = 37;

% catc = readtable('/Users/jdwulet/Google Drive/barakphase/TimeTrace10sec3z6.csv');
% second = 10;
% whichz = 6;
% z=8;
% starttime = 17;
% endtime = -1;

% catc = readtable('/Users/jdwulet/Google Drive/barakphase/TimeTrace10sec7z3.csv');
% second = 10;
% whichz = 3;
% z=8;
% starttime = -1;
% endtime = -1;

% catc = readtable('/Users/jdwulet/Google Drive/barakphase/TimeTrace10sec10z3.csv');
% second = 10;
% whichz = 3;
% z=8;
% starttime = -1;
% endtime = -1;

% catc = readtable('/Users/jdwulet/Google Drive/barakphase/TimeTrace10sec12z7.csv');
% second = 10;
% whichz = 7;
% z=8;
% starttime = -1;
% endtime = -1;

% catc = readtable('/Users/jdwulet/Google Drive/barakphase/TimeTrace10sec13z6.csv');
% second = 10;
% whichz = 6;
% z=8;
% starttime = -1;
% endtime = -1;

% catc = readtable('/Users/jdwulet/Google Drive/barakphase/TimeTrace10sec18z7.csv');
% second = 10;
% whichz = 7;
% z=8;
% starttime = 31;
% endtime = -1;

catc = readtable('/Users/jdwulet/Google Drive/barakphase/TimeTrace10sec11z7.csv');
second = 10;
whichz = 7;
z=8;
starttime = -1;
endtime = -1;


addpath('C:\Users\dwuletj\Documents\GitHub\UniversalCode\')
addpath('/Users/jdwulet/Documents/GitHub/UniversalCode/');

% R=bfopen(DataOut.Location);
% omeMeta = R{1, 4};
% stackSizeX = omeMeta.getPixelsSizeX(0).getValue(); % image width, in pixels
% stackSizeY = omeMeta.getPixelsSizeY(0).getValue(); % image height, in pixels
%
% voxelSizeX = omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER); % in µm
% voxelSizeXdouble = voxelSizeX.doubleValue(); % The numeric value represented by this object after conversion to type double
% voxelSizeY = omeMeta.getPixelsPhysicalSizeY(0).value(ome.units.UNITS.MICROMETER); % in µm
% voxelSizeYdouble = voxelSizeY.doubleValue(); % The numeric value represented by this object after conversion to type double
%
% PixelSize = voxelSizeXdouble; % in um

for i=1:size(catc,2)-2
    
    timecourse = catc(:,i);
    getname = fieldnames(timecourse);
    timecourse = timecourse.(getname{1});
    timecourse = timecourse(whichz:z:end);
   
    if starttime == -1
        st=1;
    else
        st = starttime;
    end

    if endtime == -1
        ed=length(timecourse);
    else
        ed=endtime;
    end
    
    timecouse = timecourse(st: ed);
    
    s = timecourse; %use calcium to determine phase
    t = 1:second:length(s)*second; %time vector
    Ts = mean(diff(t)); %% Sampling period
    Fs = 1/Ts;  % Sampling Frequency
    L = length(s);
    NFFT=2^nextpow2(L);
    TF=Fs./2*linspace(0,1,NFFT/2+1);
    TF(1)=[];
    TF=1./TF;
    Fn = Fs/2;  %Nyquist Frequency
    L  = length(s); %length of signal
    
    fts = fft(s-mean(s));
    
    course =s-mean(s);
    f = Fs*(0:(ceil(L/2)))/L;
    
    P2 = abs(fts/L);
    P1 = P2(1:(ceil(L/2)+1));
    P1(2:end-1) = 2*P1(2:end-1);
    
    amp_fts = abs(fts);      % Spectrum Amplitude
    phs_fts = angle(fts);
    P1ordered = sort(P1);
    highestfreq = f(find(P1==max(P1)));%this is the prob
    value = find(f==highestfreq);
    freqvec(i) = f(value);
    ampvec(i) = amp_fts(value);
    phasevec(i) = unwrap(phs_fts(value));
    
    if i==1
        FunFreq = highestfreq;
    end
    
    value1 = find(f==FunFreq);
    freqvec_Fun(i) = f(value1);
    ampvec_Fun(i) = amp_fts(value1);
    phasevec_Fun(i) = unwrap(phs_fts(value1));
    phasevec1(i) = angle(fts(value1));
    
    
    %figure
    plot(timecourse-mean(timecourse));
    hold on
%     plot(max(P1)*sin(2*pi*freqvec(i)*t));
%     hold on 
    plot(P1(value1)*sin(2*pi*freqvec_Fun(i)*t));
    hold on
    
end
phasevec1=phasevec1./(2*pi*FunFreq);
% TF(value1)
% freqvec_Fun

[minphase, minind] = min(phasevec1);
[maxphase, maxind] = max(phasevec1);
%%additional after calculate distance

difphase = maxphase-minphase;
DataOut.difphase = difphase;
DataOut.minind = minind;
DataOut.maxind = maxind;

% DataOut.maxdifphase = max(DataOut.difphase);
% DataOut.allspeed = DataOut.dist./DataOut.difphase
% DataOut.maxdist = DataOut.dist(find(DataOut.difphase == max(DataOut.difphase)));
% DataOut.speed = DataOut.allspeed(find(DataOut.difphase == DataOut.maxdifphase));
