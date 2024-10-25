function [ tail2P ] = RunFtestandTtest( x1, x2, s1, s2, n1, n2, semboolean )
%UNTITLED2 runs Ftest and if it fails the ftest it runs the ttest with
%unequal variances

%set sem
if semboolean 
    sem1 = s1;
    sem2 = s2;
    s1 = sem1*sqrt(n1);
    s2 = sem2*sqrt(n2);
end

v1 = n1-1;
v2 = n2-1;

F = s1^2/s2^2;

if s1 > s2
    p = fcdf(F,v1,v2,'upper');
else 
    p = fcdf(F,v1,v2);
end

if p > .025
    s = sqrt(((n1-1)*s1^2+(n2-1)*s2^2)/(n1+n2-2));
    s1 = s;
    s2 = s;
end

t = (x1-x2)/sqrt(s1^2/n1+s2^2/n2);

v = n1+n2-2;
tail2P = 2*tcdf(-abs(t),v);
tdist2T = @(t,v) (1-betainc(v/(v+t^2),v/2,0.5));    % 2-tailed t-distribution
T2 = [1-tdist2T(t,v)  tail2P];

end

