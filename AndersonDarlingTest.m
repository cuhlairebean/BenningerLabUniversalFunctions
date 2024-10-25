x = [0;1;0.283;0];

[h,p] = adtest(x);

y = [1;0;0;0;1;1;0;0;0;1;0;0];
%z = [x;y;1;1;0;0;0;0;0;];
    
[hy,py] = adtest(y);
%[hz,pz] = adtest(z);

p = sum(y)/length(y);
q=1-p;
stdev = sqrt(p*q);
standarddeviation = std(y);

