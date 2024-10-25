function [ Images ] = registerImages( Images )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%INPUT: Image stack (n x n x T) images
%OUTPUT: Registered image stack
fixedImage=Images(:,:,1);

for i=2:size(Images,3)
   movingImg=Images(:,:,i);
   
   [optimizer,metric]=imregconfig('monomodal');
   tform=imregtform(movingImg,fixedImage,'rigid',optimizer,metric);
   
   movingReg=imwarp(movingImg,tform,'OutputView',imref2d(size(fixedImage)));
   Images(:,:,i)=movingReg;
   
   display(['Registering Image ',num2str(i),' out of ',num2str(size(Images,3))]);
   
   
    
end



end

