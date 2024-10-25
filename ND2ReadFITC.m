function [ SepChannels,Cname,timeVec ] = ND2ReadV2(  )
%ND2READV2 Summary of this function goes here
%   Detailed explanation goes here
%Author: Matthew Westacott


%OUTPUTS: SepChannels- txn array of separated channels
            %Cname : The names of each channel 
            %timeVec : Time at each image

%Open ND2 File
R=bfopen();
%Get the number of channels
channels=R{4}.getChannelCount(0);
pics=R{1};
pics=pics(:,1);
Cname=[];
for i=0:channels-1
    %Get the channel name
Cname{i+1}=R{4}.getChannelName(0,i);
    %Get the current channel and order it
SepChannels(:,i+1)=pics((i+1):channels:length(pics));
end
timeVec=[];
timeVec(1)=0;
for j=1:size(SepChannels,1);
    %Get the time vector
    T=double(R{4}.getPlaneDeltaT(0,j));
    timeVec=[timeVec,T];
end


timeVec=timeVec';



end

