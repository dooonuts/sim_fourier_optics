% %%% Image - Black and White
% resolution = 25
% % figure(1);clf
% % hold on
% [x_unit, y_unit] = circle(0,0,1,resolution);
% % plot(x_unit,y_unit)
% fill(x_unit,y_unit,'k');
% 

paper = zeros(360,360)
[y x]= size(paper)
r = 120

for i=1:y
    for j=1:x
        if ((i-y/2)^2)+((j-x/2)^2)<(r^2);  %define origin is at the center
            paper(i,j) = 1;  %define array inside the circle eq. = 1
        end
    end
end
% figure(1);clf
% imshow(paper)

% hold off

% Low Res Sample Creation from High Res

% Illumination Plane Pattern Generation  - What does that look like at the object
x = linspace(0,2*pi,100);
y = 1+ cos(x);
size(y)
ones_axis = zeros(100);
size(ones_axis)
image = [y, ones_axis]

figure(2);clf
imshow(image)
% Multiplication of Pattern and Sample in Regular Space = Complex Wave Sure 

% Convert Complex Wave into Frequency space

% Convert PSF into Frequency (just delta in spatial domain, constant in
% frequency domain)

% Multiply the two

% Inverse Fourier Transform of image 

% Repeat all steps with phase shaped inputs