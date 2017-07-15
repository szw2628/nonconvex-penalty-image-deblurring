function [psnrvalue]=psnr(original,test)

% PSNR=10*log_10(L^2/MSE)
% MSE=1/MN*sum_i=1^Msum_j=1^n(test(i,j)-original(i,j))^2
%计算原始图像的信号功率
% A=imread(original);
%A=rgb2gray(A);
A=double(original);
% B=imread(test);
%B=rgb2gray(B);
B=double(test);

%计算MSE
%判断输入图像是否有效
[m,n]=size(A);
[m2,n2]=size(B);
if m2~=m||n2~=n
    error('图像选择错误');
end

%计算MSE
msevalue=0;
for i=1:m
    for j=1:n
        msevalue=msevalue+(A(i,j)-B(i,j))^2;
    end
end
msevalue=msevalue/(m*n);
if msevalue==0
    error('图像选择错误');
end

%计算信噪比，峰值信噪比
psnrvalue=255^2/msevalue;
psnrvalue=10*log10(psnrvalue);
end