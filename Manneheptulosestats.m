%manneheptoluse statistics
clear all
close all
clc
%% control
x1 = .718;
sem1 = 0.048948;
n1 = 24;

x2 = .606;
sem2 = 0.032897;
n2 = 24;

controlp = RunFtestandTtest( x1, x2, sem1, sem2, n1, n2, 1 );
controlp

%% 3MH
x1 = .575;
sem1 = 0.04112;
n1 = 10;

x2 = 0.211;
sem2 = 0.078368;
n2 = 10;

MH3p = RunFtestandTtest( x1, x2, sem1, sem2, n1, n2, 1 );
MH3p

%% 5MH
x1 = 0.055;
sem1 = 0.026424;
n1 = 8;

x2 = 0.217;
sem2 = 0.082803;
n2 = 8;

MH5p = RunFtestandTtest( x1, x2, sem1, sem2, n1, n2, 1 );
MH5p

%% 10MH
x1 = 0.083;
sem1 = 0.035764;
n1 = 6;

x2 = 0.061;
sem2 = 0.020239;
n2 = 6;

MH10p = RunFtestandTtest( x1, x2, sem1, sem2, n1, n2, 1 );
MH10p
