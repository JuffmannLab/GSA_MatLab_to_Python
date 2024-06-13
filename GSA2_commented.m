%% Description:
% This script performs laser wavefront shaping using a specified target image
% and laser field. It implements an iterative algorithm for hologram
% generation on a Spatial Light Modulator (SLM) and visualizes the
% resulting phase and intensity distributions.


% Clear workspace and close all figures
clear all; 
close all;
tic; % Start timing the script

% Add path to the library
addpath('.\lib');

%% Define dimensions of the image and SLM
Nx = 1920; % Number of pixels along x-axis
Ny = 1152; % Number of pixels along y-axis

%% Specify SLM (Spatial Light Modulator) dimensions
lambda = 1.035e-6; % Wavelength in meters
dx = 9.2e-6; % Pixel size in x-direction
dy = 9.2e-6; % Pixel size in y-direction
L1 = Nx * dx; % Length of x-side
L2 = Ny * dy; % Length of y-side

% Generate coordinate grid
x = -L1/2 : dx : L1/2 - dx; % Source coordinates in x
y = -L2/2 : dy : L2/2 - dy; % Source coordinates in y
[X, Y] = meshgrid(x, y);

%% Load and prepare target image
%Target = rgb2gray(imread('PSI_red_TARGET.png')); 
Target = rgb2gray(imread('2gauss.png')); % Load target image and convert to grayscale. 
% NOTE: Image size in pixel must be an even number
[p, q] = size(Target);

% Define region of interest (ROI) for the target image
shiftr = 0; % Vertical shift
shiftu = 0; % Horizontal shift
ROIX = Nx/2 - q/2 - shiftu : Nx/2 + q/2 - shiftu - 1;
ROIY = Ny/2 - p/2 - shiftr : Ny/2 + p/2 - shiftr - 1;

% Initialize target image in the ROI
Targetn = zeros(Ny, Nx);
Targetn(Ny/2 - p/2 : Ny/2 + p/2 - 1, Nx/2 - q/2 : Nx/2 + q/2 - 1) = abs(255 - Target);

%Targetn(ROIX,ROIY)=Target;
Targetn = Targetn / sum(sum(Targetn)); % Normalize target intensity
Targetn = Targetn.^(1/2); % Square root of intensity for electric field

%% Define Gaussian Laser field
x0 = 0; % Center x-coordinate
y0 = 0; % Center y-coordinate
sigma1 = 2.7e-3 / 2 * 3; % Sigma value for Gaussian
sigma2 = 2.7e-3 / 2 * 3; % Sigma value for Gaussian
Amp = 1; % Amplitude of the Gaussian

% Calculate intensity distribution
res = (2 * (X - x0).^2 / sigma1^2 + 2 * (Y - y0).^2 / sigma1^2);
Intensity = (Amp * exp(-res));
Intensity = Intensity / sum(sum(Intensity)); % Normalize intensity
Efield = Intensity.^(1/2); % Electric field amplitude

% Define mask for the target region
mask = zeros(Ny, Nx);
mask(ROIY, ROIX) = 1;

% Initialize SLM and mask constraints
SLM = ones(Ny, Nx);
Tf(1:Ny, 1:Nx) = 1;
Tf(sqrt(X.^2 + Y.^2) >= 5.3e-3) = 0; % Circular mask constraint

% Initial phase (zero) for the target
Phase_TEN = zeros(Ny, Nx);
D = Targetn .* exp(1i * Phase_TEN); % Initial deflection

% Apply mask constraint to the electric field
Efield = Efield .* Tf;

% Parameters for iterative algorithm
error = [];
iteration_num = 60; % Number of iterations
z = 0.1; % Propagation distance

%% Iterative optimization using Gerchberg-Saxton Algorithm (GSA)
for i = 1:iteration_num
    tic;
    % Propagate backward and display
    A = propagateback(D, -1, lambda, dx, dy, Nx, Ny, X, Y, z, z, z);
    A = A ./ sqrt(sum(sum(abs(A).^2))); % Normalize field amplitude
    as = subplot(3,5,1); imagesc(abs(A).^2); title('SLM ideal Int.');
    as = [as subplot(3,5,6)]; imagesc(angle(A)); title('SLM ideal Phase');

    % Apply amplitude constraint to the hologram plane and display
    B = Efield .* exp(1i * angle(A));
    as = [as subplot(3,5,2)]; imagesc(abs(B).^2); title('SLM true Int.');
    as = [as subplot(3,5,7)]; imagesc(angle(B)); title('SLM true Phase');

    % Propagate forward and display
    %C = 1/sqrt(Nx*Ny)*fftshift(fft2(fftshift(B)));
    C = propagateforward(B, 1, lambda, dx, dy, Nx, Ny, X, Y, z, z, z);
    C = C ./ sqrt(sum(sum(abs(C).^2))); % Normalize field amplitude
    as = [as subplot(3,5,3)]; imagesc(abs(C).^2); title('Image Int.');
    as = [as subplot(3,5,8)]; imagesc(angle(C)); title('Image Phase');

    % Update and display deflected field with target phase
    D = Targetn .* exp(1i * angle(C)); %noise window
    as = [as subplot(3,5,4)]; imagesc(abs(D).^2); title('Image ideal Int.');
    as = [as subplot(3,5,9)]; imagesc(angle(D)); title('Image enforced phase');

    % Compute and display error and efficiency
    error = sum(sum(abs(abs(C .* mask).^2 - abs(Targetn .* mask).^2)));
    deff = sum(sum(abs(C .* mask).^2)) / sum(sum(abs(C).^2));
    subplot(3,5,5); hold on; plot(i, error, 'bo'); hold off; title('Error');
    subplot(3,5,10); hold on; plot(i, deff, 'bo'); hold off; title('Efficiency');

    toc;
    drawnow;
    linkaxes(as);
end

% Final visualizations
figure(3211); imagesc(abs(C).^2); % Display final intensity
figure(1111); imagesc(angle(C)); % Display final phase

%% Generate image for display on the SLM
data = angle(B); % Hologram phase data

% Normalize the data to the range [0, 1]
min_val = min(data(:));
max_val = max(data(:));
normalized_data = (data - min_val) / (max_val - min_val);

% Convert normalized data to 8-bit unsigned integers
image_data = uint8(normalized_data * 255);

% Save the image as BMP
imwrite(image_data, '2gauss-marius.bmp');
