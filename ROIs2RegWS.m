function [sRegions] = ROIs2RegWS(cvsROIs, vnImageSize)

% Modified from ROIs2Regions by D. Muir
% Modified by W. Schleicher - April 2020
% ROIs2RegWS - FUNCTION Convert a set of imported ImageJ ROIs into a Matlab regions structure
%
% Usage: [sRegions] = ROIs2Regions(cvsROIs, vnImageSize, numcells)
%
% 'cvsROIs' is a cell array of ImageJ ROI structures, as imported by
% ReadImageJROI. 'vnImageSize' is a vector [M N] containing the size of the
% image in pixels. 'numcells' is number of ROIs in the archive.  
%
% 'sRegions' will be a structure compatible with the Matlab regions
% structure format, as returned by bwconncomp. It will contain one region
% for each compatible ROI in 'cvsROIs'.
%
% Only a subset of ImageJ ROI types is supported for conversion:
% 'rectangle', oval', 'polygon' and 'freehand'.

% Original Author: Dylan Muir <muir@hifo.uzh.ch>
% Created: 2011
% Modified by Wolfgang Schleicher <wolfgang.schleicher@cuanschutz.edu>
% Modified April 2020

% - Check arguments

if (nargin < 2)
   disp('*** ROIs2Regions: Incorrect usage.');
   help ROIs2Regions;
   return;
end

% - Build a regions structure
sRegions.Connectivity = 8;
sRegions.ImageSize = vnImageSize;
sRegions.NumObjects = numel(cvsROIs);
CellMask = zeros(vnImageSize);
% sRegions.PixelIdxList = {};

for nROIIndex = 1:numel(cvsROIs)
   sThisROI = cvsROIs{nROIIndex};
   
   switch (lower(sThisROI.strType))
      case 'rectangle'
         if (isfield(sThisROI, 'strSubtype') && isequal(lower(sThisROI.strSubtype), 'shape'))
            % - Skip this one
            continue;
            
         else
            % - Make a rectangular mask
                ROIMask = false(vnImageSize);
                sThisROI.vnRectBounds = sThisROI.vnRectBounds + 1;
                ROIMask(sThisROI.vnRectBounds(1):sThisROI.vnRectBounds(3), sThisROI.vnRectBounds(2):sThisROI.vnRectBounds(4)) = true;
%             sRegions.PixelIdxList{nROIIndex} = find(mbThisMask');
         end
         
      case 'oval'
         % - Draw an oval inside the bounding box
             ROIMask = ellipse2mask('bounds', vnImageSize, sThisROI.vnRectBounds+1);
             CellMask = CellMask + ROIMask.*nROIIndex;
             CellMask(find(CellMask>nROIIndex)) = nROIIndex;
%          mbThisMask = ellipse2mask('bounds', vnImageSize, sThisROI.vnRectBounds+1);
%          sRegions.PixelIdxList{nROIIndex} = find(mbThisMask');
         
      case {'polygon'; 'freehand'}
         % - Draw a polygonal mask
             ROIMask = poly2mask(sThisROI.mnCoordinates(:,1)+1,sThisROI.mnCoordinates(:,2)+1,vnImageSize(1),vnImageSize(2));
             CellMask = CellMask + ROIMask.*nROIIndex;
             CellMask(find(CellMask>nROIIndex)) = nROIIndex;
%          mbThisMask = poly2mask(sThisROI.mnCoordinates(:, 1)+1, sThisROI.mnCoordinates(:, 2)+1, vnImageSize(1), vnImageSize(2));
%          sRegions.PixelIdxList{nROIIndex} = find(mbThisMask');
         
      otherwise
         warning( 'ROIs2Regions:unsupported', ...
                  '--- ROIs2Regions: Warning: Unsupported ROI type.');
   end
end
sRegions.CellMask = CellMask;

% --- END of ROIs2RegWS.m ---
