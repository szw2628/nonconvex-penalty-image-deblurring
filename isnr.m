function [isnrvalue]=isnr(original,Blured_img,test)

% ISNR=10*log_10(||blured_img-original||^2/||test-original||^2)


A=double(original);

B=double(test);

%º∆À„MSE
%≈–∂œ ‰»ÎÕºœÒ «∑Ò”––ß
[m,n]=size(A);
[m2,n2]=size(B);
if m2~=m||n2~=n
    error('ÕºœÒ—°‘Ò¥ÌŒÛ');
end

%º∆À„MSE
msevalue=0;
for i=1:m
    for j=1:n
        msevalue=msevalue+(A(i,j)-B(i,j))^2;
    end
end
%msevalue=msevalue/(m*n);
if msevalue==0
    error('ÕºœÒ—°‘Ò¥ÌŒÛ');
end

%solve ISNR
signal=0;
for i=1:m
    for j=1:n
        signal=signal+(A(i,j)-Blured_img(i,j))^2;
    end
end
%signal=signal/(m*n);
isnrvalue=signal/msevalue;
isnrvalue=10*log10(isnrvalue);
end