%% Description
% This script performs a laser wavefront shaping algorithm using Gaussian
% approximation and other transformations on optical fields. It includes 
% propagation calculations, creation of masks, and application of Zernike 
% polynomials for aberration correction. The script also generates visualizations 
% of the transformed fields and evaluates the results.

% Clear workspace and close all figures
clear all; 
close all;
tic; % Start timing the script

% Add path to the library
addpath('.\lib');

%% Define dimensions of the input grid (image dimentions)
Nx = 2446; % Number of pixels along x-axis
Ny = 2446; % Number of pixels along y-axis

%% Define the region of interest (ROI) for the target image

Nimg = 60; % Number of pixels in the target image.
global shiftr shiftu
shiftr = 0; % Vertical shift
shiftu = 0; % Horizontal shift

ROIX = Nx/2 - Nimg/2 - shiftu : Nx/2 + Nimg/2 - shiftu;
ROIY = Ny/2 - Nimg/2 - shiftr : Ny/2 + Nimg/2 - shiftr;

%% Define mask for SLM (Spatial Light Modulator)
Nimg_x = 1152;
Nimg_y = 1920;
ROI_spx = Nx/2 - Nimg_x/2 - shiftu : Nx/2 + Nimg_x/2 - shiftu - 1;
ROI_spy = Ny/2 - Nimg_y/2 - shiftr : Ny/2 + Nimg_y/2 - shiftr - 1;
mask_sp = zeros(Ny, Nx);
mask_sp(ROI_spx, ROI_spy) = 1;

ROIX2 = Nx/2 - 1148/2 - 1 : Nx/2 + 1148/2;
ROIY2 = Ny/2 - 1148/2 - 1 : Ny/2 + 1148/2;

mask = zeros(Ny, Nx);
mask(ROIY, ROIX) = 1;

mask2 = zeros(Ny, Nx);
mask2(ROIY2, ROIX2) = 1;

lambda = 1.035e-6; % Wavelength in meters

% Define spatial resolution
dx = 9.2e-6; % Pixel size in x-direction
dy = 9.2e-6; % Pixel size in y-direction
L1 = Nx * dx; % Length of x-side
L2 = Ny * dy; % Length of y-side

x = -L1/2 : dx : L1/2 - dx; % Source coordinates in x
y = -L2/2 : dy : L2/2 - dy; % Source coordinates in y
[X, Y] = meshgrid(x, y);
x0 = 0;
y0 = 0;

%% Define Gaussian distribution at the SLM
constant = 1;
sigma1 = constant * 2.7e-3 / 2 * 3;
sigma2 = constant * 2.7e-3 / 2 * 3;

Amp = 1;
res = 2 * ((X - x0).^2 / sigma1^2 + (Y - y0).^2 / sigma1^2);
Intensity = Amp * exp(-res);
Intensity = Intensity / sum(sum(Intensity)); % Normalize intensity
Efield = Intensity.^(1/2); % Electric field amplitude

%% Define target Gaussian distribution
constant = 0.00800; % Target Gaussian constant. Pondemorotive deflection Nimg=60
%  constant=0.08; %dummy deflection Nimg=600
sigma1 = constant * 2.7e-3;
res = 2 * ((X - x0).^2 / sigma1^2 + (Y - y0).^2 / sigma1^2);
Intensity = Amp * exp(-res);
Intensity = Intensity / sum(sum(Intensity)); % Normalize intensity
EfieldTarget = Intensity.^(1/2); % Target electric field amplitude
figure(321321);
imagesc(EfieldTarget); % Display target electric field

%% Define noise region and initialize target
Target = rgb2gray(imread('PSI_red_TARGET.JPG'));
Target = imresize(Target, [Ny Nx]);
Target = double(Target);
%Target=zeros*(Target); % commented 2024
NX = 1920;
NY = 1152;

% TARGET = imresize(Target, [NY NX]); % commented 2024, alternative to the
% next row
TARGET = imresize(Target*0, [NY NX]); % Initialize TARGET to zero
%TARGET = zeros(NY, NX); % Other way to initialize TARGET to zero
Zernikepadded = TARGET; % Initialize Zernike pattern arrays
Zernikepadded_tilt = TARGET;
Zernikepadded_defocus = TARGET;
Zernikepadded_spherical = TARGET;
Zernikepadded_astx = TARGET;
Zernikepadded_comay = TARGET;
Zernikepadded_asty = TARGET;
Zernikepadded_comax = TARGET;
Zernikepadded_helical = TARGET;

%% Define signal region using an image
apple = imread('SinusoidGrating_Period128.bmp');
apple = imrotate(apple, 270);
apple = imresize(apple, [Nimg + 1, Nimg + 1]);
apple = double(apple);
apple = apple / sum(sum(apple)); % Normalize image
MAX = max(apple(:));

%% Initialize SLM and mask constraints
SLM = ones(Nx, Ny);
[nor, noc] = size(SLM);
Tf(1:nor, 1:noc) = 1;
maskconstraint = Tf;
maskconstraint_gaussian = Tf;
Tf(sqrt(X.^2 + Y.^2) >= 2.3e-3) = 0;
Tf2(1:nor, 1:noc) = 0;
Tf2(sqrt(X.^2 + Y.^2) >= 5.3e-3) = 1;
maskconstraint(sqrt((X + shiftu * dx).^2 + (Y + shiftr * dx).^2) >= Nimg/2 * dx + 0.2e-5) = 0; % Wavefront shaping mask

maskconstraint_gaussian(sqrt((X + shiftu * dx).^2 + (Y + shiftr * dx).^2) >= Nimg/10 * dx + 0.2e-5) = 0;
%Target(mask==1)=EfieldTarget(Ny/2-Nimg/2:Ny/2+Nimg/2,Nx/2-Nimg/2:Nx/2+Nimg/2);%gaussian; commented 2024
%Target=Target.*maskconstraint_gaussian; % commented 2024

const = 0;
D = Target .* exp(1i * const); % Initial deflection

Efield = Efield .* Tf;

%% Iterative optimization using Gerchberg-Saxton Algorithm (GSA)
k = 1;
K = 1;
error = [];
iteration_num = 1;
clims = [0 MAX];

Z_p = 0.5; % Sampling parameter. Note: non critical sampling - <0.2
for i = 1:iteration_num
    tic;
    %A = sqrt(Nx*Ny)*ifftshift(ifft2(ifftshift(D)));
    A = propagateback1(D, -1, lambda, dx, dy, Nx, Ny, X, Y, Z_p, Z_p, Z_p);
    A = A ./ sqrt(sum(sum(abs(A).^2))) .* Tf;
%       as = subplot(3,5,1); imagesc(abs(A).^2); title('SLM ideal Int.'); 
%   	as = [as subplot(3,5,6)]; imagesc(angle(A)); title('SLM ideal Phase'); 

    B = Efield .* exp(1i * angle(A)) .* Tf; % Apply amplitude constraint to the hologram plane
    %sum(sum(abs(B).^2))
%   	as = [as subplot(3,5,2)]; imagesc(abs(B).^2); title('SLM true Int.'); 
%   	as = [as subplot(3,5,7)]; imagesc(angle(B)); title('SLM true Phase'); 

    %C = 1/sqrt(Nx*Ny)*fftshift(fft2(fftshift(B)));
    C = propagateforward1(B, 1, lambda, dx, dy, Nx, Ny, X, Y, Z_p, Z_p, Z_p);
    C = C ./ sqrt(sum(sum(abs(C).^2))) .* Tf;
%   	as = [as subplot(3,5,3)]; imagesc(abs(C).^2); title('Image Int.'); 
%   	as = [as subplot(3,5,8)]; imagesc(angle(C)); title('Image Phase'); 

    D = C; % Update D with propagated field
    D = K * abs(Target) .* exp(1i * angle(C)) .* Tf; % Apply GSA
    %D=D/sqrt(sum(sum(abs(D).^2)));
%   	as = [as subplot(3,5,4)]; imagesc(abs(D).^2); title('Image ideal Int.');
%   	as = [as subplot(3,5,9)]; imagesc(angle(D)); title('Image enforced phase');

    
    % Compute error and correlation
    error = sum(sum(abs(abs(C).^2 - abs(Target).^2))); %60 px image
    correlationCODE = corr2(abs(C .* mask).^2, abs(Target .* mask).^2);
    deff = sum(sum(abs(C .* mask).^2)) / sum(sum(abs(C).^2));
%   	subplot(3,5,5); hold on ; plot(i,error,'bo'); hold off; title('Error');
%   	subplot(3,5,10); hold on ; plot(i,deff,'bo'); hold off; title('Efficiency');
%   	subplot(3,5,11); hold on ; plot(i,correlation,'bo'); hold off; title('Correlation');

    
    toc;
%   drawnow;
%	linkaxes(as);
%   pause;
end

figure(5);
imagesc(abs(C(1215:1235, 1214:1234)).^2);
%figure();subplot(1,2,1); imagesc(abs(C).^2); subplot(1,2,2); imagesc(abs(Target).^2);


%% Define Zernikes in appropriate dimensions
ROIX2 = NX/2 - 1148/2 - 1 : NX/2 + 1148/2;
ROIY2 = NY/2 - 1148/2 - 1 : NY/2 + 1148/2;
mask4 = zeros(NY, NX);
mask4(ROIY2, ROIX2) = 1;

% Create and apply Zernike polynomials
k = 5; 
DX = 0.00174;
Zernike = CreateZernikePolynomials(DX, k);
Zernikepadded_defocus(mask4 == 1) = Zernike;
Zernikepadded_defocus(isnan(Zernikepadded_defocus)) = 0;
amp_defocus = 0.1:0.1:1;

k = 9; 
DX = 0.00174;
Zernike = CreateZernikePolynomials(DX, k);
Zernikepadded_tilt(mask4 == 1) = Zernike;
Zernikepadded_tilt(isnan(Zernikepadded_tilt)) = 0;
amp_tilt = 0.1:0.1:1;

k = 13; 
DX = 0.00174;
Zernike = CreateZernikePolynomials(DX, k);
Zernikepadded_spherical(mask4 == 1) = Zernike;
Zernikepadded_spherical(isnan(Zernikepadded_spherical)) = 0;
amp_sph = 0:0.1:2;

k = 4; 
DX = 0.00174;
Zernike = CreateZernikePolynomials(DX, k);
Zernikepadded_astx(mask4 == 1) = Zernike;
Zernikepadded_astx(isnan(Zernikepadded_astx)) = 0;

k = 6; 
DX = 0.00174;
Zernike = CreateZernikePolynomials(DX, k);
Zernikepadded_asty(mask4 == 1) = Zernike;
Zernikepadded_asty(isnan(Zernikepadded_asty)) = 0;

k = 8; 
DX = 0.00174;
Zernike = CreateZernikePolynomials(DX, k);
Zernikepadded_comax(mask4 == 1) = Zernike;
Zernikepadded_comax(isnan(Zernikepadded_comax)) = 0;

k = 9; 
DX = 0.00174;
Zernike = CreateZernikePolynomials(DX, k);
Zernikepadded_comay(mask4 == 1) = Zernike;
Zernikepadded_comay(isnan(Zernikepadded_comay)) = 0;

%% Propagate Gaussian with spherical aberrations
Z_sp = zeros(Ny, Nx);
Z_sp(mask_sp == 1) = Zernikepadded_spherical;
% t1 = B.*exp(1i*(1*Z_sp));      
t1 = exp(1i * (-1 * Z_sp));
t2 = propagateforward1(t1, 1, lambda, dx, dy, Nx, Ny, X, Y, 0.2, 0.2, 0.2);
t3 = propagateforward1(B, 1, lambda, dx, dy, Nx, Ny, X, Y, 0.2, 0.2, 0.2);
figure(2121); 
subplot(1, 2, 1);
imagesc(abs(t2(1215:1235, 1214:1234)).^2);
subplot(1, 2, 2);
imagesc(abs(t3(1215:1235, 1214:1234)).^2);

%% Propagate helical pattern
angl = angle(A);
Nxs = 2446; Nys = Nxs;
dxs = 1;
dys = dxs;
L1s = Nxs * dxs; % Side length
L2s = Nys * dys;
xs = -L1s/2 : dxs : L1s/2 - dxs;
ys = -L2s/2 : dys : L2s/2 - dys;
[Xs, Ys] = meshgrid(xs, ys);
helicalpattern = atan2(Ys, Xs) .* Tf;
% Zernikepadded_helical(mask4==1) = helicalpattern; 
coef = 3;
B2 = Efield .* exp(1i * (angl)) .* exp(1i * (coef * helicalpattern));
ab2 = angle(B2);
heli = exp(1i * (coef * helicalpattern));
ah = angle(heli);
C2 = propagateforward1(B2, 1, lambda, dx, dy, Nx, Ny, X, Y, Z_p, Z_p, Z_p);
f = 0.2;
x1 = linspace(-lambda * f / 2 / dx, lambda * f / 2 / dx, Nx);
y1 = linspace(-lambda * f / 2 / dy, lambda * f / 2 / dy, Ny);
figure(71);
subplot(2, 3, 1);
imagesc(angl(860:1600, 860:1600));
title('Hologram with bb');
subplot(2, 3, 2);
imagesc(ah(860:1600, 860:1600));
title('Vortex 3');
subplot(2, 3, 3);
imagesc(ab2(860:1600, 860:1600));
title('Superposition h V');
subplot(2, 3, 4);
imagesc(abs(C2(1160:1290, 1160:1290)).^2);
title('Intensity');
subplot(2, 3, 5);
imagesc(abs(C2(1160:1290, 1160:1290)).^2);
title('Fit');

% figure(444); imagesc(abs(C2(1160:1290,1160:1290)).^2);axis equal; colormap(bluewhitered);
% print(gcf,'x4dep.png','-dpng','-r300');

C3 = abs(C2(1160:1290, 1160:1290)).^2;
line = C3(65, :);

%% Fit polynomial to line data
xf = linspace(1, length(line), length(line));
yf = line;
% y = y/136*2*pi;

ind = xf >= 49 & xf <= 81;

p = polyfit(xf(ind), yf(ind), 4); 
f = polyval(p, xf(ind)); 
plot(xf, yf, 'o', xf(ind), f, '-'); 
legend('Data', 'x^4');
