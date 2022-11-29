% %%% Image - Black and White
% resolution = 25
% % figure(1);clf
% % hold on
% [x_unit, y_unit] = circle(0,0,1,resolution);
% % plot(x_unit,y_unit)
% fill(x_unit,y_unit,'k');
% 

pixels = 360;

sample_image = zeros(pixels,pixels);
[y x]= size(sample_image);
r = 120

for i=1:y
    for j=1:x
        if ((i-y/2)^2)+((j-x/2)^2)<(r^2);  %define origin is at the center
            sample_image(i,j) = 1;  %define array inside the circle eq. = 1
        end
    end
end
% figure(1);clf
% imshow(paper)

% hold off


%% adding random code to test
phase
% Low Res Sample Creation from High Res

% Illumination Plane Pattern Generation  - What does that look like at the object
lambda = 500e-3;
NA = 0.95;
k0 = 2*NA/lambda;
cos_wave = cos_image_generation(30,k0);
% cos_wave = ones(pixels,pixels);
% imshow(cos_wave);

% Multiplication of Pattern and Sample in Regular Space = Complex Wave 
complex_sample = cos_wave;
% complex_sample = cos_wave.*sample_image;
complex_sample2 = cos_wave2.*sample_image;
% imshow(complex_sample)

% Convert Complex Wave into Frequency space
freq_complex_sample = fft2(complex_sample);
% freq_complex_sample2 = fft2(complex_sample2);
% m = abs(freq_complex_sample)
% phase = unwrap()
% imshow(freq_complex_sample)
% figure(3);clf
% reconstructed_image = ifft2(freq_complex_sample);
% imshow(reconstructed_image)

% Multiply with PSF in Frequency Domain (Incoherent OTF for Circular Aperture)

w = pixels;
z2 = 1;

fx = linspace(0,w-1,w);
fy = linspace(0,w-1,w);
[FX,FY] = meshgrid(fx,fy);

rho = FX.^2 + FY.^2;
rho_0 = w/2/lambda/z2;
H_incoh_freq = 2/pi * (acos(rho/2/rho_0) - rho/2/rho_0 * sqrt(1 - (rho/2/rho_0)^2));

E_fft = I1_fft .* H_incoh_freq;

% Inverse Fourier Transform of image 

x_lin = linspace(0,10,pixels);
y_lin = cos(x_lin*pi/4)
y_fft = fft(y_lin)
xdft = fftshift(y_fft);
plot(xdft,y_fft)

%% Repeat all steps with phase shaped inputs
% for i=1:3
%     phase = 60*(i-1);
%     cos_wave = cos_image_generation(phase)
%     complex_sample = cos_wave.*sample_image;
%     freq_complex_sample = fft2(complex_sample)
% %     k0 = 
% 
% 
% end


