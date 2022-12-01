close all
clear
clc

%% Ground truth of sample

% Circle for now
pixels = 512;
sample = zeros(pixels,pixels); %create empty array
[y,x] = size(sample); %define y,x as size of array
r = 100; %define radius of a circle
for i=1:y
    for j=1:x
        if ((i-y/2)^2)+((j-x/2)^2)<(r^2)  %define origin is at the center
            sample(i,j) = 1;  %define array inside the circle eq. = 1
        end
    end
end

% figure; imshow(sample)
%% Create sample (low resolution of sample)

%% Generate illumination excitation pattern

lambda = 500e-9;
% lambda = 75;
NA = 0.95;
k0 = 2*NA/lambda;
% k0 = 75;
theta_ = [0 60 120];
phi_ = [0 120 240];
phi_ = deg2rad(phi_);

pattern_pixels = pixels*2;
x = 1:pattern_pixels; 
y = 1:pattern_pixels; 
[X,Y] = meshgrid(x,y);

for t = 1:length(theta_)
        theta = theta_(t);
    for p = 1:length(phi_)
        illum = 1+cos(2*pi*k0/pixels*X+phi_(p));
        rotated_illum = imrotate(illum,theta);
        win = centerCropWindow2d(size(rotated_illum),size(sample));
        rotated_illum = imcrop(rotated_illum,win);
%         
%         figure; 
%         imshow(rotated_illum);

    %% Multiply sample and pattern in real space, initial intensity I1
        I1 = rotated_illum.*sample;
%         I1 = sample;
%         figure; imshow(I1)

        I1_(:,:,t) = I1;
        %% Convert I1 to frequency space

        I1_fft = fftshift(fft2(I1));
%         I1_fft = fft2(I1);
%         figure; imagesc(log(1+abs(I1_fft)))
%         title("FFT after pattern")

        I1_fft_(:,:,t) = I1_fft;
        %% Multiply with PSF (OTF)
        w = pixels;

        fx = linspace(0,w-1,w);
        fy = linspace(0,w-1,w);
        [FX,FY] = meshgrid(fx,fy);

    %         z2 = 1;
    %     rho = FX.^2 + FY.^2;
    %     rho_0 = w/2/lambda/z2;
    %     H_incoh_freq = 2/pi * (acos(rho/2/rho_0) - rho/2/rho_0 * sqrt(1 - (rho/2/rho_0)^2));
    %     H_incoh_freq = sample;

        scale = 0.63;
        R=sqrt(min(FX,abs(FX-w)).^2+min(FY,abs(FY-w)).^2);
        yy=abs(2*besselj(1,scale*R+eps,1)./(scale*R+eps)).^2;
        OTF2d=fft2(yy);
        OTF2dmax = max(max(abs(OTF2d)));
        OTF2d = OTF2d./OTF2dmax;
        OTF2dc = abs(fftshift(OTF2d));
%         OTF2dc = abs(OTF2d);
        H_incoh_freq = OTF2dc;

        E_fft = I1_fft .* H_incoh_freq;
%         figure; imagesc(log(1+abs(E_fft)))
% %         figure; imshow(E_fft)
%         title("FFT after OTF")

%         E_fft_(:,:,t) = E_fft;
        %% Inverse FT of output intensity, I2

        E = ifft2(ifftshift(E_fft));
%         figure; imagesc(E)
%         title("IFFT after OTF")

        E_(:,:,t,p) = E_fft;
    end
end

M = [
    [1 1/2*exp(-1i*deg2rad(0)) 1/2*exp(1i*deg2rad(0))]
    [1 1/2*exp(-1i*deg2rad(120)) 1/2*exp(1i*deg2rad(120))]
    [1 1/2*exp(-1i*deg2rad(240)) 1/2*exp(1i*deg2rad(240))]
];
% t1 = [E_(:,:,1,1) E_(:,:,1,2) E_(:,:,1,3)];
Minv = inv(M);

ft1o = Minv(1,1)*E_(:,:,1,1) + Minv(1,2)*E_(:,:,1,2) + Minv(1,3)*E_(:,:,1,3);
ft1p = Minv(2,1)*E_(:,:,1,1) + Minv(2,2)*E_(:,:,1,2) + Minv(2,3)*E_(:,:,1,3);
ft1m = Minv(3,1)*E_(:,:,1,1) + Minv(3,2)*E_(:,:,1,2) + Minv(3,3)*E_(:,:,1,3);

ft2o = Minv(1,1)*E_(:,:,2,1) + Minv(1,2)*E_(:,:,2,2) + Minv(1,3)*E_(:,:,3,3);
ft2p = Minv(2,1)*E_(:,:,2,1) + Minv(2,2)*E_(:,:,2,2) + Minv(2,3)*E_(:,:,3,3);
ft2m = Minv(3,1)*E_(:,:,2,1) + Minv(3,2)*E_(:,:,2,2) + Minv(3,3)*E_(:,:,3,3);


%% Shift Frequency components by theta

% 1,: = 0
% 2,: = 60
% 3,: = 120

% theta, phi
% 0, 0
% 0, 120 
% 0, 240


%% doubling Fourier domain size if necessary
w=pixels;
% DoubleMatSize = 0;
% if ( 2*Kotf > wo )
% 	DoubleMatSize = 1; % 1 for doubling fourier domain size, 0 for keeping it unchanged
% end
% if ( DoubleMatSize>0 )
%     t = 2*w;
%     to = t/2;
%     u = linspace(0,t-1,t);
%     v = linspace(0,t-1,t);
%     [U,V] = meshgrid(u,v);
%     fDoTemp = zeros(2*w,2*w);
%     fDpTemp = zeros(2*w,2*w);
%     fDmTemp = zeros(2*w,2*w);
%     OTFtemp = zeros(2*w,2*w);
%     fDoTemp(wo+1:w+wo,wo+1:w+wo) = fDof;
%     fDpTemp(wo+1:w+wo,wo+1:w+wo) = fDpf;
%     fDmTemp(wo+1:w+wo,wo+1:w+wo) = fDmf;
%     OTFtemp(wo+1:w+wo,wo+1:w+wo) = OTFo;
%     clear fDof fDpf fDmf OTFo w wo x y X Y
%     fDof = fDoTemp;
%     fDpf = fDpTemp;
%     fDmf = fDmTemp;
%     OTFo = OTFtemp;
%     clear fDoTemp fDpTemp fDmTemp OTFtemp
% else
%%
t = w;
to = t/2;
u = linspace(0,t-1,t);
v = linspace(0,t-1,t);
[U,V] = meshgrid(u,v);

Est_1 = fft2(ifft2(E_(:,:,1,1)));
k1 = k0/t.*(U-to) + k0/t.*(V-to);
k2 = k0*cosd(60)/t.*(U-to) + k0*1i*sind(60)/t.*(V-to);
k3 = k0*cosd(120)/t.*(U-to) + k0*1i*sind(120)/t.*(V-to);
% k0_ref = fft2(ifft2(E_(:,:,2,1)).*exp(+1i.*2*pi*k0));
% figure;imshow(k0_ref)

Est_2 = fft2(ifft2(ft1p).*exp(+1i.*2*pi*k2));
% Est_3 = fft2(ifft2(E_(:,:,3,1)).*exp(+1i.*2*pi*k3));
% Est_4 = fft2(ifft2(E_(:,:,1,1)));
% Est_5 = fft2(ifft2(E_(:,:,2,1)).*exp(+1i.*2*pi*(120/t.*(U-to) + 60/t.*(V-to))));
% Est_6 = fft2(ifft2(E_(:,:,3,1)).*exp(-1i.*2*pi*(120/t.*(U-to) + 120/t.*(V-to))));
% Est_4 = fft2(ifft2(E_))

Est_sum = Est_2;
figure; imagesc(log(1+abs(ft1p)));
% figure;imagesc(log(1+abs(Est_sum)));
% figure;imshow(Est_sum)
% Est_sum = Est_6 + Est_5;
% figure;imshow(Est_sum)
% Est_1 = fft2( ifft2(E_(:,:,1,1)).*exp(-1j*0) );
% Est_2 = fft2( ifft2(E_(:,:,2,1)).*exp(-1j*deg2rad(60)) );
% Est_3 = fft2( ifft2(E_(:,:,3,1)).*exp(-1j*deg2rad(120)) );
% 
% Est_4 = fft2( ifft2(E_(:,:,1,2)).*exp(+1j*0) );
% Est_5 = fft2( ifft2(E_(:,:,2,2)).*exp(+1j*deg2rad(60)) );
% Est_6 = fft2( ifft2(E_(:,:,3,2)).*exp(+1j*deg2rad(120)) );


%% Merge using Weiner Filter

% w = size(pixels,1);
% wo = w/2;
% fx = linspace(0,w-1,w);
% fy = linspace(0,w-1,w);
% 
% % [fX,fY] = meshgrid(fx,fy);
% 
% S1 = conj(H_incoh_freq+k0)*Est_1 / abs(H_incoh_freq+k0)^2;
% S2 = conj(H_incoh_freq+k0)*Est_2 / abs(H_incoh_freq+k0)^2;
% S3 = conj(H_incoh_freq+k0)*Est_3 / abs(H_incoh_freq+k0)^2;
% 
% S4 = conj(H_incoh_freq+k0)*Est_4 / abs(H_incoh_freq+k0)^2;
% S5 = conj(H_incoh_freq+k0)*Est_5 / abs(H_incoh_freq+k0)^2;
% S6 = conj(H_incoh_freq+k0)*Est_6 / abs(H_incoh_freq+k0)^2;
% 
% S7 = conj(H_incoh_freq+k0)*Est_7 / abs(H_incoh_freq+k0)^2;
% S8 = conj(H_incoh_freq+k0)*Est_8 / abs(H_incoh_freq+k0)^2;
% S9 = conj(H_incoh_freq+k0)*Est_9 / abs(H_incoh_freq+k0)^2;
% 
% 
% % S = S1+S2+S3+S4+S5+S6+S7+S8+S9;
% 
% S = Est_1 + Est_2 + Est_3 + Est_4 + Est_5 + Est_6 + Est_7 + Est_8 + Est_9;


% figure; imshow(S);
% title("Reconstruction in Frequency Domain")
% %%
% output = ifft2(S);
% 
% figure; imshow(output);
% title("IFFT of reconstruction")
