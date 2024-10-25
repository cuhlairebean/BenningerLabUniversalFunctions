function [ amplitudeVal ] = HighGAnalyze( data )
%HIGHGANALYZE Summary of this function goes here
%   Detailed explanation goes here
amplitudeVal=[];

   dMap=data.drivingMap;
   
   amplitudeVal=[amplitudeVal;nanmean(dMap(find(dMap)))];
    
    



end

