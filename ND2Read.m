% Thomas Hraha <thomas.hraha@ucdenver.edu> or <tomhraha@gmail.com>
%
% function to read proprietary image formats into a 3D matrix in matlab
% it calls: bfUpgradeCheck.m, bfsave.m, bfOpen3DVolume.m, bfopen.m,
% bfGetReader.m, bfGetPlan.m, bfGetFileExtensions.m, bfCheckJavaPath.m,
% imreadBF.m, imreadBFmeta.m and loci_tools.jar either directly or
% indirectly. 
%
% the function will call a window where you will need to select the FOLDER 
% which contains the image file.
% 
% input note: filename must include file extension (e.g., .nd2) and be 
% input as a string. for example:
%
% [stack] = ND2Read('14April_PI_w200_2mm_r3.nd2')
%
%%

function [stack] = ND2Read(filename)

% select directory
P = get(findobj('Tag','DirectoryList'),'UserData');
    if ~isdir(P)
        P = fileparts(mfilename('fullpath'));
    end
        
% make sure it is valid
    currentFile = uigetdir(P,'Select Folder');
    if currentFile == 0
        return
        elseif ~isdir(currentFile)
            return;
    end
currentFile = ([currentFile '\' filename]);

meta=imreadBFmeta(currentFile); clc;
[stack] = imreadBF(currentFile,meta.zsize,1:meta.nframes,meta.channels); 
clc;
     
end
