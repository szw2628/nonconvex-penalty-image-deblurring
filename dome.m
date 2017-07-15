%%%%%%    min_x  1/2 || Ax - y ||_2^2 + lambda * phi( Wx , a) %%%%%%%
%%%%%%               + gamma * phi( Dx , b)                   %%%%%%%
clear
clc
x = double(imread('barbara256.tif'));
%x = imresize(x,[256 256]);
[m, n] = size(x);
imshow(uint8(x));
%%%%%    create a gaussian blurring operator   %%%%%
mask_type = 'gaussian';
mask_size = 11;
mask_deviation = 1.6;
mask = fspecial(mask_type, mask_size, mask_deviation);
%%%%%    create a blurred image with a GWN     %%%%%
 %standard_deviation = sqrt(2);
 standard_deviation =2;
Blurred_img = imfilter(x, mask, 'circular');
Blurred_img = Blurred_img + standard_deviation*randn(m, n);
figure;
imshow(uint8(Blurred_img));
%%%%%          some known parameters           %%%%%
%  lambda  comes from  lambda * phi( Wx , a)   %
lambda = 0.0095;%lambda=0.0095
%  gamma   comes from  gamma * phi( Dx , b)    %
gamma = 0.55;%gamma=0.0055
%  mu      comes from  mu/2 * ||u - v + d||_2^2    %
mu = 0.0001;%mu=0.0001
%  rho     comes from  rho/2 * ||u - u^||_2^2   %
%          and rho/2 * ||alpha - alpha^||_2^2   %
rho = 0.0001;%rho=0.0001
%  eta     comes from  eta/2 * ||Wu - alpha + d_1||_2^2  %
eta = 0.0003;%eta=0.002
%  tau     comes from   tau/2 * ||D_{h}v - d_x +b_x||_2^2  %
%          and tau/2 * ||D_{v}v - d_y +b_y||_2^2 %
tau = 0.003;%tau=0.002
%    a     comes from lambda * phi( Wx , a)    %
a = 0.99*(eta + rho)/lambda;
%    b     comes from lambda * phi( Dx , b)    %
b = 0.99*tau/gamma;
%    l   is length of framelet
l = 5;
%    Nit is times of interation
Nit = 30;

%%%%%compute eigen_value diagonal matrix of mask%%%%%
Sbig = psf2otf(mask, [m, n]);
ATA = conj(Sbig).*Sbig;
ATy = conj(Sbig).*fft2(Blurred_img);
%%%%compute eigen_value diagonal matrixof gradient%%%
D = [-1, 1];
DBig = psf2otf(D, [m, n]);
DDBig = psf2otf(D', [m, n]);
DhTDh=conj(DBig).*DBig;
DvTDv=conj(DDBig).*DDBig;
%%%%%%%%%   some inintional  conditions    %%%%%%%%%%
u_new = Blurred_img;% u0=blurred image%
v_new = Blurred_img;% v0=blurred image%
d_new = 0 * Blurred_img;% d0=0，d为Lagrange常数%
alpha_new = QASDC(Blurred_img, 1, 2);%alpha0是小波分解的结果，小波系数%
d = QASDC(0 * Blurred_img, l, 2);%Wu=d新的变量%
b_x = 0 * Blurred_img;%初始值b_x=0%
b_y = 0 * Blurred_img;%初始值b_y=0%
U = u_new;
f=zeros(1,Nit);
%%%%%% algorithm for solve a stable point   %%%%%%%%%
tic;
kk=0;
%kk=1;
%     z = Sbig .* x - fft2(Blurred_img);
%     %ztz = conj(z).* z;
%     %f1 = 1/2 * trace(ztz);
%     f1 = 1/2 * norm(z,2) * norm(z,2);
%     
%     wu = QASDC(x ,1 ,2);
%     f2 = 0;
%     for i = 1 : l
%         for j = 1 : l
%             phi{i,j} =  2/(a*sqrt(3)) * (atan((abs(wu{i,j}).*a*2+1)/sqrt(3)) - pi/6);
%             f2 = f2 + lambda * norm(phi{i,j},1);
%         end
%     end
%     %phiwu = QASRE(phi ,1 ,2);
%     %f2 = lambda * norm(phiwu ,1);
%     
%     DhX = imfilter(x,D,'circular');
%     DvX = imfilter(x,D','circular');
%     Dx = sqrt((DhX)^2 +(DvX)^2);
%     phidx = 2/(a*sqrt(3)) * (atan((abs(Dx).*a*2+1)/sqrt(3)) - pi/6);
%     f3 = gamma * norm(phidx ,1);
%     
%     f(:,kk) = f1 + f2 + f3;
for k = 1 : Nit
    kk=kk+1;
    u_old = u_new;
    v_old = v_new;
    d_old = d_new;
    alpha_old = alpha_new;
    %%%%%   Mean Doubly Augmemted Lagrangian(MDAL) method   %%%%%
    for i = 1:l
        for j = 1:l
            temp{i, j} = alpha_old{i,j} - d{i,j};
        end
    end
    TEMP = QASRE(temp, 1, 2);
    u_new = ifft2((ATy + fft2(mu * (v_old - d_old) + eta * TEMP + rho * u_old))./(ATA + mu + eta + rho));
    u_new = real(u_new);
    Temp = QASDC(u_new, 1, 2);
    for i = 1:l
        for j = 1:l
            TEmp = (eta * (Temp{i, j} + d{i, j}) + rho * (alpha_old{i, j}))/(eta + rho);
            alpha_new{i, j} = thresh_atan(TEmp, lambda/(eta + rho), a);
        end
    end
%     alpha_new{1, 1} = (eta * (Temp{1, 1} + d{1, 1}) + rho * (alpha_old{1, 1}))/(eta + rho);
    alpha_new{1, 1} = Temp{1, 1};
    for i = 1:l
        for j = 1:l
            d{i, j} = d{i, j} + 1 * (Temp{i, j} - alpha_new{i, j});
        end
    end
    U = (k + 1)/(k + 2) * U + 1/(k + 2) * u_new;
    %%%%%   The split Bregman algorithm of isotropic proximal operator    %%%%%
    d_x = imfilter(v_old,D,'circular');
    d_y = imfilter(v_old,D','circular');
    v_new = ifft2((mu * fft2(u_new + d_old) + tau * (conj(DBig) .* fft2(d_x - b_x) +conj(DDBig) .* fft2(d_y - b_y)))./(mu + tau * (DhTDh + DvTDv)));
    v_new = real(v_new);
    Dhv = imfilter(v_new,D,'circular');
    Dvv = imfilter(v_new,D','circular');
    s = sqrt((Dhv + b_x).^2 + (Dvv + b_y).^2);
    d_x = thresh_atan(s, gamma/tau, b);
    d_x = d_x .* ((Dhv + b_x)./s);
    d_y = thresh_atan(s, gamma/tau, b);
    d_y = d_y .* ((Dvv + b_y)./s);
    b_x = b_x + 1 * (Dhv - d_x);
    b_y = b_y + 1 * (Dvv - d_y);
    %%%%%   the three subproblem   %%%%%
    d_new = d_old + 0.01 * (u_new -v_new);
    
    %%%%    the energy of objective function  %%%%%
    z = Sbig .* U - fft2(Blurred_img);
    %ztz = conj(z).* z;
    %f1 = 1/2 * trace(ztz);
    f1 = 1/2 * norm(z,2) * norm(z,2);
    
    wu = QASDC(U ,1 ,2);
    f2 = 0;
    for i = 1 : l
        for j = 1 : l
            phi{i,j} =  2/(a*sqrt(3)) * (atan((abs(wu{i,j}).*a*2+1)/sqrt(3)) - pi/6);
            f2 = f2 + lambda * norm(phi{i,j},1);
        end
    end
    %phiwu = QASRE(phi ,1 ,2);
    %f2 = lambda * norm(phiwu ,1);
    
    DhU = imfilter(U,D,'circular');
    DvU = imfilter(U,D','circular');
    Du = sqrt((DhU)^2 +(DvU)^2);
    phidu = 2/(a*sqrt(3)) * (atan((abs(Du).*a*2+1)/sqrt(3)) - pi/6);
    f3 = gamma * norm(phidu ,1);
    
    f(:,kk) = f1 + f2 + f3;
   
end

    m = 0 : 1 : 29;
    plot (m,f(1,:));
toc;
%%%%  test %%%%
Deblurred_img = real(U);
Deblurred(Deblurred_img > 255) = 255; Deblurred_img(Deblurred_img < 0) = 0;
psnrvalue1=psnr(x,Blurred_img);
psnrvalue2=psnr(x,Deblurred_img);
isnrvalue=isnr(x,Blurred_img,Deblurred_img);
figure;
imshow(uint8(Deblurred_img));
% imshow(uint8(Deblurred_img),'Border','tight');