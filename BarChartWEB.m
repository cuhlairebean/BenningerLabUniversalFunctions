% Tom Hraha <Thomas.Hraha@ucdenver.edu> OR <Tomhraha@gmail.com>
% 
% Adapted from:
% 'http://www.mathworks.com/support/solutions/en/data/1-18H2E/?solution=1-18H2E'
%
%
%%

function [varargout] = BarChartWEB(model_series, model_error)
% Code for Figure 1 period robustness
%a = [1.5 2.86; 3 3.144];
%b = [0.245 0.41; 0.6264 0.6831];
a = model_series;
b = model_error;
gw = 0.75 ; % almost touching, 1 is touching
x = 1:size(a,1) ;
figure
fig(1) = bar(x-gw/4,a(:,1),gw/2,'w') ;
hold on ;
fig(2) = bar(x+gw/4,a(:,2),gw/2,'b') ;
h1 = errorbar(x-gw/4,a(:,1),b(:,1),'k','LineStyle','none');
h2 = errorbar(x+gw/4,a(:,2),b(:,1),'k','LineStyle','none');
set(fig(1),'FaceColor',[1 1 1]) ;
set(fig(2),'FaceColor',[0 0 0]) ;


hb = get(h1,'children'); 
Xdata = get(hb(2),'Xdata');
temp = 4:3:length(Xdata);
temp(3:3:end) = [];
% xleft and xright contain the indices of the left and right
% endpoints of the horizontal lines
xleft = temp; xright = temp+1; 
% Increase line length by 0.2 units
Xdata(xleft) = Xdata(xleft) - .05;
Xdata(xright) = Xdata(xright) + .05;
set(hb(2),'Xdata',Xdata)

hb = get(h2,'children'); 
Xdata = get(hb(2),'Xdata');
temp = 4:3:length(Xdata);
temp(3:3:end) = [];
% xleft and xright contain the indices of the left and right
% endpoints of the horizontal lines
xleft = temp; xright = temp+1; 
% Increase line length by 0.2 units
Xdata(xleft) = Xdata(xleft) - .05;
Xdata(xright) = Xdata(xright) + .05;
set(hb(2),'Xdata',Xdata)
end