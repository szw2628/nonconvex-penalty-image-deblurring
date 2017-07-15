function imag = QASRE(coeffs, L,types)
%Reconstruction 1 level by quai-affine system
%%%%%%%%%%%%20-07-2013 by Shen ZhengWei%%%%%%%%%%%%%%%%%%%%%%%
%   Detailed explanation goes here
if types==1
    %reconstruction filters which is the time-inverse of decomposition
    %filter
    
     h{1}=(1/4)*[1,2,1];
    h{2}=(sqrt(2)/4)*[-1,0,1];
    h{3}=(1/4)*[-1,2,-1];
   
elseif types==2
    h{1}=(1/16)*[1,4,6,4,1];
    h{2}=(1/8)*[1,2,0,-2,-1];
    h{3}=(sqrt(6)/16)*[1,0,-2,0,1];
    h{4}=(1/8)*[1,-2,0,2,-1];
    h{5}=(1/16)*[1,-4,6,-4,1];
end



ImSize=size(coeffs{1,1});
R=zeros(ImSize);
NoF=length(h);
imag=zeros(ImSize);

% [mc,nc]=size(coeffs);
for i=1:L
    
    for ir=1:NoF
       
        for jr=1:NoF
             temp=zeros(ImSize);
             %first columun convlution
             temp=imfilter(coeffs{ir,jr}',h{1,jr},'circular');
             %second row convlution
             R=R+imfilter(temp', h{1,ir},'circular'); 
        end 
             
    end
end
imag=R;

end

