function [ cMat,count ] = ConnectionMat( CT,count )
%CONNECTIONMAT Summary of this function goes here
%   Detailed explanation goes here
for k=1:size(CT,1)
numCoord(k)=length(find(CT(k,:)));
end

SCT=size(CT,1);
cMat=zeros(SCT,SCT);

[y,x]=max(numCoord);
if ~isempty(x)
useMax=find(CT(x,:));
cMat(x,CT(x,useMax))=count;
vals=CT(x,useMax);
cMat(x,x)=count;
for i=1:length(useMax)
cMat(vals(i),vals(i))=count;
end
count=count+1;
end

[C,G]=sort(numCoord);
G=fliplr(G);
for i=2:size(CT,1)
    i2=G(i);
    if i2~=x && cMat(i2,i2)==0
        cMat(i2,i2)=count;
   usePos=CT(i2,find(CT(i2,:))); 
    if isempty(usePos)
        cMat(i2,i2)=count;
        count=count+1;
    end
    if ~isempty(usePos)
       for j=1:length(usePos)
          if  cMat(usePos(j),usePos(j))==0
             cMat(i2,2)=count;
             cMat(i2,usePos(j))=count; 
             count=count+1;
          else
             
             % cMat(i2,i2)=cMat(usePos(j),i2);
              
          end
       end

        
    end
   
   
   
end



end
count=count+1;

end

