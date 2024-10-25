function [DataOut]=ImgPrep(img)

matrixClusterLo=img;
[MSize,NSize,TSize]=size(matrixClusterLo); 
tt=1:TSize;
 
 %% Correcting for negative intensity values
 
%              disp('checking for negative intensity values');
%             
%             
%             act_min=zeros(1,TSize);act_max=zeros(1,TSize);
% 
%             
%             for k=1:TSize
%                 act_min(k)=(min(min(matrixClusterLo(:,:,k),[],2),[],1));
%                 act_max(k)=(max(max(matrixClusterLo(:,:,k),[],2),[],1));
%                 
%                 %check for negative intensity values
%                 
%                 for i=1:MSize
%                     for j=1:NSize
%                         
%                         correct the negative values
%                         if matrixClusterLo(i,j,k) < 0
%                             disp('negative values');
%                             
%                             matrixClusterLo(i,j,k) = act_max(k) + (abs(act_min(k)) + obj.imgStack(i,j,k));
%                             
%                         else
%                             obj.imgStack(i,j,k) = images(i,j,k);
%                             
%                         end
%                         
%                     end
%                 end
%             end
% 


 
 %% Averaging filter applied to every stack
 
 h = fspecial('average', 10);
%  for i=1:TSize
% avgIslet(:,:,i)=imfilter(matrixClusterLo(:,:,i),h,'replicate'); 
%  end
 avgIslet=imfilter(matrixClusterLo,h,'replicate');
%  size(avgIslet)
% figure
% imshow(avgIslet(:,:,1))
% close
 DataOut=avgIslet;
    %% Correcting for background intensity
%  
%  disp('select background area');
% 
%  
%  figure
% I=mat2gray(double(mean(avgIslet(:,:,:),3)));
% imshow(I,[min(min(I)) max(max(I))*(.75)])
% 
%  disp('correcting for background intensity');
%  
% bkgMask=createMask(imfreehand);  %ones(MSize,NSize);  
% bkgMaskStack = repmat(bkgMask,[1 1 TSize]);
% bkgmatrixCluster = avgIslet.* bkgMaskStack;
% 
% bkgThresh=mean(shiftdim(mean(mean(bkgmatrixCluster,2),1),1))
% 
% for m=1:MSize
%     for n=1:NSize
%         
%         
%         if avgIslet(m,n) <= bkgThresh
%             
%             avgIslet(m,n,:)=0;
%             
% %             avgIsletOld(m,n)=avgIslet(m,n);
%             
%             for t=1:TSize
%             avgIslet(m,n,t)=0;
% 
%             end
%             
% %             avgIsletNew(m,n)=avgIslet(m,n);
%             
%         end
%             
%       
%         
%     end
% end
% 
% DataOut=avgIslet;

%  DataOut.avgIsletOld=avgIsletOld;
%  DataOut.avgIsletNew=avgIsletNew;
%% Correcting for photobleaching
        
% mn_avgIslet=mean(avgIslet,3);
%  for m=1:MSize
%      for n=1:NSize
% 
% 
% 
%          if avgIslet(m,n)>0
%              
%          course=shiftdim(smooth(shiftdim(avgIslet(m,n,start:endd),1)),1);
%          dataFit  = polyfit(tt,course,2);
%          dataFit2 = polyval(dataFit,tt);
%          
%          
%         courseNorm=(course./dataFit2)-ones(length(tt),1)';
%         
%         matrixNorm(m,n,1:length(tt))=courseNorm;
%         
%         
%         matrixStd(m,n)=std(courseNorm);
%         
%          end
%  
%      end
%  end
%  
%  avgstdIslet=std(matrixNorm,0,3)./mean(mean(std(matrixNorm,0,3),2),1);
%  
%   avgmatrixNorm=matrixNorm; %.*img;
%   
%    DataOut.StdAvgIslet=matrixStd;
%  DataOut.AvgMatrixNorm=avgmatrixNorm;
%  
  %%
  

 %%
 
 
