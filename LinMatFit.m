function [ A ] = LinMatFit( A )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

T=1:size(A,1);

for i=1:size(A,2)
X=A(:,i);

[FIT]=polyfit(T,X',1);

f=polyval(FIT,T);

A(:,i)=A(:,i)-f';
    
end

end

