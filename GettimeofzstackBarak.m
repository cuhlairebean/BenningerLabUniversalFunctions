clear all
close all
clc
T = readtable('/Users/jdwulet/Google Drive/barakphase/Timecourses/Time11.csv');
cell1 = T.cell1;
cell2 = T.cell2;
cell3 = T.cell3;
cell4 = T.cell4;
realtime = T.realtime;
z=[1,3,6,8];
zstacks = 8;

cell1 = cell1(z(1):zstacks:end);
cell2 = cell2(z(2):zstacks:end);
cell3 = cell3(z(3):zstacks:end);
cell4 = cell4(z(4):zstacks:end);

cell1 = cell1/mean(cell1);
cell2 = cell2/mean(cell2);
cell3 = cell3/mean(cell3);
cell4 = cell4/mean(cell4);

realtime = realtime(z(1):zstacks:end);

cells=[cell1,cell2,cell3,cell4];

%% 15
T = readtable('/Users/jdwulet/Google Drive/barakphase/Timecourses/Time15.csv');
cell1 = T.cell1;
cell2 = T.cell2;
cell3 = T.cell3;
cell4 = T.cell4;
z=[1,2,4,8];
zstacks = 8;

cell1 = cell1(z(1):zstacks:end);
cell2 = cell2(z(2):zstacks:end);
cell3 = cell3(z(3):zstacks:end);
cell4 = cell4(z(4):zstacks:end);

cell1 = cell1/mean(cell1);
cell2 = cell2/mean(cell2);
cell3 = cell3/mean(cell3);
cell4 = cell4/mean(cell4);

cells=[cell1,cell2,cell3,cell4];

%% 7
T = readtable('/Users/jdwulet/Google Drive/barakphase/Timecourses/Time7.csv');
cell1 = T.cell1;
cell2 = T.cell2;
cell3 = T.cell3;
cell4 = T.cell4;
z=[2,2,5,8];
zstacks = 8;

cell1 = cell1(z(1):zstacks:end);
cell2 = cell2(z(2):zstacks:end);
cell3 = cell3(z(3):zstacks:end);
cell4 = cell4(z(4):zstacks:end);

cell1 = cell1/mean(cell1);
cell2 = cell2/mean(cell2);
cell3 = cell3/mean(cell3);
cell4 = cell4/mean(cell4);

cells=[cell1,cell2,cell3,cell4];

%% 45
T = readtable('/Users/jdwulet/Google Drive/barakphase/Timecourses/Time45.csv');
cell1 = T.cell1;
cell2 = T.cell2;
cell3 = T.cell3;
cell4 = T.cell4;
z=[2,7,7,1];
zstacks = 12;

cell1 = cell1(z(1):zstacks:end);
cell2 = cell2(z(2):zstacks:end);
cell3 = cell3(z(3):zstacks:end);
cell4 = cell4(z(4):zstacks:end);

cell1 = cell1/mean(cell1);
cell2 = cell2/mean(cell2);
cell3 = cell3/mean(cell3);
cell4 = cell4/mean(cell4);

cells=[cell1,cell2,cell3,cell4];

%% 22
T = readtable('/Users/jdwulet/Google Drive/barakphase/Timecourses/Time22.csv');
realtime = T.realTime;

cell1 = T.cell1;
cell2 = T.cell2;
cell3 = T.cell3;
cell4 = T.cell4;
z=[1,1,2,3];
zstacks = 3;

realtime = realtime(z(1):zstacks:end);

cell1 = cell1(z(1):zstacks:end);
cell2 = cell2(z(2):zstacks:end);
cell3 = cell3(z(3):zstacks:end);
cell4 = cell4(z(4):zstacks:end);

cell1 = cell1/mean(cell1);
cell2 = cell2/mean(cell2);
cell3 = cell3/mean(cell3);
cell4 = cell4/mean(cell4);

cells=[cell1,cell2,cell3,cell4];

%% 33
T = readtable('/Users/jdwulet/Google Drive/barakphase/Timecourses/Time33.csv');
%realtime = T.realTime;

cell1 = T.cell1;
cell2 = T.cell2;
cell3 = T.cell3;
cell4 = T.cell4;
z=[1,1,1,1];
zstacks = 1;

%realtime = realtime(z(1):zstacks:end);

cell1 = cell1(z(1):zstacks:end);
cell2 = cell2(z(2):zstacks:end);
cell3 = cell3(z(3):zstacks:end);
cell4 = cell4(z(4):zstacks:end);

cell1 = cell1/mean(cell1);
cell2 = cell2/mean(cell2);
cell3 = cell3/mean(cell3);
cell4 = cell4/mean(cell4);

cells=[cell1,cell2,cell3,cell4];

%% 35
T = readtable('/Users/jdwulet/Google Drive/barakphase/Timecourses/Time35v3.csv');
%realtime = T.realTime;

cell1 = T.cell1;
cell2 = T.cell2;
cell3 = T.cell3;
%cell4 = T.cell4;
z=[1,1,1,1];
zstacks = 1;

%realtime = realtime(z(1):zstacks:end);

cell1 = cell1(z(1):zstacks:end);
cell2 = cell2(z(2):zstacks:end);
cell3 = cell3(z(3):zstacks:end);
%cell4 = cell4(z(4):zstacks:end);

cell1 = cell1/mean(cell1);
cell2 = cell2/mean(cell2);
cell3 = cell3/mean(cell3);
%cell4 = cell4/mean(cell4);

cells=[cell1,cell2,cell3];%,cell4];



%% 35
T = readtable('/Users/jdwulet/Google Drive/barakphase/Timecourses/Time38.csv');
%realtime = T.realTime;

cell1 = T.cell1;
cell2 = T.cell2;
cell3 = T.cell3;
cell4 = T.cell4;
z=[1,1,1,1];
zstacks = 1;

%realtime = realtime(z(1):zstacks:end);

cell1 = cell1(z(1):zstacks:end);
cell2 = cell2(z(2):zstacks:end);
cell3 = cell3(z(3):zstacks:end);
cell4 = cell4(z(4):zstacks:end);

cell1 = cell1/mean(cell1);
cell2 = cell2/mean(cell2);
cell3 = cell3/mean(cell3);
cell4 = cell4/mean(cell4);

cells=[cell1,cell2,cell3,cell4];

%% 4
T = readtable('/Users/jdwulet/Google Drive/barakphase/Timecourses/Time42.csv');
%realtime = T.realTime;

cell1 = T.cell1;
cell2 = T.cell2;
cell3 = T.cell3;
cell4 = T.cell4;
cell5 = T.cell5;
z=[1,4,4,7,4];
zstacks = 8;

%realtime = realtime(z(1):zstacks:end);

cell1 = cell1(z(1):zstacks:end);
cell2 = cell2(z(2):zstacks:end);
cell3 = cell3(z(3):zstacks:end);
cell4 = cell4(z(4):zstacks:end);
cell5 = cell5(z(5):zstacks:end);

cell1 = cell1/mean(cell1);
cell2 = cell2/mean(cell2);
cell3 = cell3/mean(cell3);
cell4 = cell4/mean(cell4);
cell5 = cell5/mean(cell5);

cells=[cell1,cell2,cell3,cell4,cell5];

%% 61127
T = readtable('/Users/jdwulet/Google Drive/barakphase/Timecourses/Time61127.csv');
%realtime = T.realTime;

cell1 = T.cell1;
cell2 = T.cell2;
cell3 = T.cell3;
cell4 = T.cell4;
z=[2,4,5,8];
zstacks = 8;

%realtime = realtime(z(1):zstacks:end);

cell1 = cell1(z(1):zstacks:end);
cell2 = cell2(z(2):zstacks:end);
cell3 = cell3(z(3):zstacks:end);
cell4 = cell4(z(4):zstacks:end);

cell1 = cell1/mean(cell1);
cell2 = cell2/mean(cell2);
cell3 = cell3/mean(cell3);
cell4 = cell4/mean(cell4);

cells=[cell1,cell2,cell3,cell4];


