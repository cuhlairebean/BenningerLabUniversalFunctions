function [ fura,chr2,timeVec ] = ND2ReadV2( filename )
%ND2READV2 Summary of this function goes here
%   Detailed explanation goes here

R=bfopen(filename);
channels=R{4}.getChannelCount(0);
pics=R{1};
pics=pics(:,1);

for i=0:channels-1
   C1=R{4}.getChannelName(0,i);
   if C1=='Fura 380'
       furaChannel=i+1;
   end
   if C1=='Fitc'
       chr2Channel=i+1;
   end
end

fura=pics(furaChannel:channels:length(pics));
chr2=pics(chr2Channel:channels:length(pics));
timeVec=[];
for j=0:channels:length(pics)-1
    T=double(R{4}.getPlaneDeltaT(0,j));
    timeVec=[timeVec,T];
end




end

