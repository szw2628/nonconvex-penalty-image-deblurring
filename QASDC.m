function outimg = QASDC(img,L,types)
%img--input image

%types----'1' is haar; '2' is B_{2}(x)with 3 filter; '4' is cubic 
%L ---level of decomposition
%UNTITLED Summary of this function goes here
%   Detailed explanation goes 

%%%%%%%%%%%%20-07-2013 by Shen ZhengWei%%%%%%%%%%%%%%%%%%%%%%%
[m,n]=size(img);
if types==1
    %dc filters
    
    h{1}=(1/4)*[1,2,1];
    h{2}=(sqrt(2)/4)*[1,0,-1];
    h{3}=(1/4)*[-1,2,-1];
   
elseif types==2
    h{1}=(1/16)*[1,4,6,4,1];
    h{2}=(1/8)*[-1,-2,0,2,1];
    h{3}=(sqrt(6)/16)*[1,0,-2,0,1];
    h{4}=(1/8)*[-1,2,0,-2,1];
    h{5}=(1/16)*[1,-4,6,-4,1];
end
% number of filters
NoF=length(h);
% decompose the img to L level
temp=zeros(size(img));
    


    for ir=1:NoF
        temp=imfilter(img,h{1,ir},'circular'); 
        for jr=1:NoF
          D{ir,jr}=imfilter(temp', h{1,jr}, 'circular')';         
        end
    end
    


outimg=D;

end

