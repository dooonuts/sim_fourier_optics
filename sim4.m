close all
clear
clc

%% Constants
lambda = 100e-3;
NA = 0.95;
k0 = 2*NA/lambda;
rho0 = 60;
k0 = 8*rho0;

%% Ground truth of sample

% Circle for now
pixels = 1024;
% sample = zeros(pixels,pixels); %create empty array
% [y,x] = size(sample); %define y,x as size of array
% r = 100; %define radius of a circle
% for i=1:y
%     for j=1:x
%         if ((i-y/2)^2)+((j-x/2)^2)<(r^2)  %define origin is at the center
%             sample(i,j) = 1;  %define array inside the circle eq. = 1
%         end
%     end
% end
sample = imread("USAF-1951.svg.png");
[r c] = size(sample)
sample = double(sample(r/2-pixels/2:r/2+pixels/2-1,r/2-pixels/2:r/2+pixels/2-1));
figure; imshow(sample)

%% Incoherent Transfer Function
circle = zeros(pixels,pixels); %create empty array
[y,x] = size(sample); %define y,x as size of array
for i=1:y
    for j=1:x
        if ((i-y/2)^2)+((j-x/2)^2)<(rho0^2)  %define origin is at the center
            circle(i,j) = 1;  %define array inside the circle eq. = 1
        end
    end
end
H_incoh_freq = conv2(circle,circle,'same');

%% Create diffraction limited sample (low resolution of sample)
sample_freq = fftshift(fft2(sample));
diffraction_limited_image = ifft2(ifftshift(sample_freq.*H_incoh_freq));
figure;imagesc(abs(diffraction_limited_image))
colormap gray


%% Generate illumination excitation pattern

% lambda = 500e-9;
theta_ = [0 60 120];
phi_ = [0 120 240];
phi_ = deg2rad(phi_);

pattern_pixels = pixels;
x = 1:pattern_pixels; 
y = 1:pattern_pixels; 
[X,Y] = meshgrid(x,y);

for t = 1:length(theta_)
    theta = theta_(t);

    kx = k0/pixels*cosd(theta);
    ky = k0/pixels*sind(theta);
   

    for p = 1:length(phi_)
        phi = phi_(p);
        rotated_illum = 1+cos(kx*X+ky*Y+phi);   
        ronchi_grating_freq = fftshift(fft2(rotated_illum))./(circle+1e-5);
        ronchi_grating_real = ifft2(ifftshift(ronchi_grating_freq));
%         figure; imshow(rotated_illum);
        figure; imagesc(log(1+abs(ronchi_grating_real)));


    %% Multiply sample and pattern in real space, initial intensity I1
        I1 = rotated_illum.*sample;
%         figure; imshow(I1)

        I1_(:,:,t) = I1;
        %% Convert I1 to frequency space

        I1_fft = fftshift(fft2(I1));
%         figure; imagesc(log(1+abs(I1_fft)))
%         title("FFT after pattern")

        I1_fft_(:,:,t) = I1_fft;
        %% Multiply with PSF (OTF)
        w = pixels;

        fx = linspace(0,w-1,w);
        fy = linspace(0,w-1,w);
        [FX,FY] = meshgrid(fx,fy);

        E_fft = I1_fft .* H_incoh_freq;
%         figure; imagesc(log(1+abs(E_fft)))
%         title("FFT after OTF")

%         E_fft_(:,:,t) = E_fft;
        %% Inverse FT of output intensity, I2

        E = ifft2(ifftshift(E_fft));
%         figure; imagesc(E)
%         title("IFFT after OTF")

        E_(:,:,t,p) = E_fft;
    end
end

%% Separate components using phase

M = [ 1 -exp(-1j*phi_(1)) -exp(1j*phi_(1))
    1 -exp(-1j*phi_(2)) -exp(1j*phi_(2))
    1 -exp(-1j*phi_(3)) -exp(1j*phi_(3))];

Minv = inv(M);

for i = 1:length(theta_)
    Y_tAp1 = E_(:,:,i,1);
    Y_tAp1 = Y_tAp1(:)'; % turn 2D matrix into row vector
    Y_tAp2 = E_(:,:,i,2);
    Y_tAp2 = Y_tAp2(:)';
    Y_tAp3 = E_(:,:,i,3);
    Y_tAp3 = Y_tAp3(:)';
    D = Minv* [Y_tAp1; Y_tAp2; Y_tAp3];
    D_tAp1 = reshape(D(1,:),[pixels,pixels]);
    D_tAp2 = reshape(D(2,:),[pixels,pixels]);
    D_tAp3 = reshape(D(3,:),[pixels,pixels]);
    
    D_tA(:,:,1) = D_tAp1;
    D_tA(:,:,2) = D_tAp2;
    D_tA(:,:,3) = D_tAp3;
    
    % D_tA = cat([D_tAp1] [D_tAp2] [D_tAp3]);
    
%     S_tA = D_tA;
    
    kx = k0/pixels*cosd(theta);
    ky = k0/pixels*sind(theta);
    
    % kshift = [-round(kx), round(ky)];
    % kshift = [round(kx*pixels-pixels/2), round(ky*pixels-pixels/2)];
    
    for j=1:length(phi_)
        maximum = max(max(D_tA(:,:,j))); [ind1,ind2] = find(D_tA(:,:,j) == maximum);
        current_shift = [round(ind1-pixels/2-1), round(ind2-pixels/2-1)];
        kshift(:,:,length(theta_)*(i-1)+j) = current_shift;
        S_tA(:,:,length(theta_)*(i-1)+j) = circshift(D_tA(:,:,j),-current_shift);
    end    
    % maximum = max(max(D_tAp2)); [ind1,ind2] = find(D_tAp2 == maximum);
    % kshift = [round(ind1-pixels/2), round(ind2-pixels/2)];
    % S_tAp2 = circshift(D_tAp2,-kshift);
    
    % figure; imagesc(abs(1+log(D_tAp2));
    
    % figure; imagesc(abs(1+log(S_tA(:,:,1)))); 
    % figure; imagesc(abs(1+log(S_tA(:,:,2)))); 
    % figure; imagesc(abs(1+log(S_tA(:,:,3))));
end  

S_tA_sum = zeros(pixels,pixels);

for i=1:length(theta_)*length(phi_)
    S_tA_sum = S_tA_sum + S_tA(:,:,i);
end    


%% Weiner Filter is Hconj/(|H|^2+Pnoise/Psignal), in our case... we have no additive noise so just becomes H^-1
% 
H_incoh_freq_sum = zeros(pixels,pixels);
for j=1:length(phi_)
    for i=1:length(theta_)
%        shifted_H_incoh_freq = H_incoh_freq;
        shifted_H_incoh_freq = circshift(H_incoh_freq, -kshift(:,:,length(theta_)*(i-1)+j));
        %tempG = conj(shifted_H_incoh_freq)./(abs(shifted_H_incoh_freq).^2).*S_tA(:,:,3*(i-1)+j);
        tempG = inv(shifted_H_incoh_freq).*(S_tA(:,:,length(theta)*(i-1)+j));
        G(:,:,length(theta_)*(i-1)+j) = tempG;
%         figure; imagesc(log(1+abs(S_tA(:,:,3*(i-1)+j))))
%         figure; imagesc(log(1+abs(H_incoh_freq)));
%         figure; imagesc(log(1+abs(inv(shifted_H_incoh_freq))));
%        figure; imagesc((tempG));
         H_incoh_freq_sum = H_incoh_freq_sum + shifted_H_incoh_freq;
    end
end

reconstructed_freq = S_tA_sum./(H_incoh_freq_sum + 1e-4);
% H_incoh_freq_sum_biased = H_incoh_freq_sum + ;
figure; imagesc(log(1+abs(reconstructed_freq)));

sim_image = ifft2(ifftshift(reconstructed_freq));
% G_sum = zeros(512,512);
% 
% for i=1:9
%     G_sum = G_sum + G(:,:,i);
% end
% G_sum = G_sum/9;

figure;imagesc(abs(sim_image));
colormap gray
% reconstruct = ifft(ifftshift(G_sum));
% figure; imagesc(log(1+abs(reconstruct)))
