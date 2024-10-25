clear all 
close all
clc

addpath('C:\Users\dwuletj\Documents\GitHub\UniversalCode\')

%% WT NADH
Location = 'C:\Users\dwuletj\Google Drive\LabStuff\Projects\Nurins Exper\NADHCalciumRepImages\WT\R0680_20G_I4.czi';
samefile = 0;
GFPsecond = 0;
imagenum = 6;
Scalebar_length_microm = 10; % in pixels
x_location = 200; % denotes x location in pixels
y_location = 400; % denotes y location in pixels
GetScaleBar(Location,'R0680_20G_I4',samefile, GFPsecond, imagenum, x_location, y_location, Scalebar_length_microm, 1);

%% x2/10 NADH
Location = 'C:\Users\dwuletj\Google Drive\LabStuff\Projects\Nurins Exper\NADHCalciumRepImages\x2_10\R0681_20G_I5.czi';
samefile = 0;
GFPsecond = 0;
imagenum = 6;
Scalebar_length_microm = 10; % in pixels
x_location = 110; % denotes x location in pixels
y_location = 370; % denotes y location in pixels
GetScaleBar(Location,'R0681_20G_I5',samefile, GFPsecond, imagenum, x_location, y_location, Scalebar_length_microm, 1);

%% 1x/10 NADH
Location = 'C:\Users\dwuletj\Google Drive\LabStuff\Projects\Nurins Exper\NADHCalciumRepImages\x1\R0718_20G_I4.lsm';
samefile = 0;
GFPsecond = 0;
imagenum = 6;
Scalebar_length_microm = 10; % in pixels
x_location = 150; % denotes x location in pixels
y_location = 410; % denotes y location in pixels
GetScaleBar(Location,'R0718_20G_I4',samefile, GFPsecond, imagenum, x_location, y_location, Scalebar_length_microm, 1);

%% WT calcium
Location = 'C:\Users\dwuletj\Google Drive\LabStuff\Projects\Nurins Exper\NADHCalciumRepImages\WT\R2142_11G_I3.czi';
samefile = 1;
GFPsecond = 2;
imagenum = 9;
Scalebar_length_microm = 10; % in pixels
x_location = 170; % denotes x location in pixels
y_location = 380; % denotes y location in pixels
GetScaleBar(Location,'R2142_11G_I3',samefile, GFPsecond, imagenum, x_location, y_location, Scalebar_length_microm, 2);


%% x2/10 calcium
%Location = 'C:\Users\dwuletj\Google Drive\LabStuff\Projects\Nurins Exper\NADHCalciumRepImages\x2_10\R2056_11G_I5.czi';
Location = 'C:\Users\dwuletj\Google Drive\LabStuff\Projects\Nurins Exper\NADHCalciumRepImages\x2_10\R0681_11G_I4.czi';
samefile = 0;
GFPsecond = 0;
imagenum = 140; %could do image 1
Scalebar_length_microm = 10; % in pixels
x_location = 150; % denotes x location in pixels
y_location = 420; % denotes y location in pixels
GetScaleBar(Location,'R0681_11G_I4',samefile, GFPsecond, imagenum, x_location, y_location, Scalebar_length_microm, 2);

%% x5/10 calcium
Location = 'C:\Users\dwuletj\Google Drive\LabStuff\Projects\Nurins Exper\NADHCalciumRepImages\x1\R0686_11G_I5.czi';
samefile = 1;
GFPsecond = 2;
imagenum = 1;
Scalebar_length_microm = 10; % in pixels
x_location = 150; % denotes x location in pixels
y_location = 390; % denotes y location in pixels
GetScaleBar(Location,'R0686_11G_I5',samefile, GFPsecond, imagenum, x_location, y_location, Scalebar_length_microm, 2);

