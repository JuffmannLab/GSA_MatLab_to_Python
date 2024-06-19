%% MATLAB Script for Generating and Modifying Holographic Patterns
% This script generates and manipulates holographic patterns, including
% setting up target regions, applying Gaussian masks, and simulating phase constraints.



%% Structure
% 1. Define dimensions of image and SLM
% 2. Define ROI for target image
% 3. 
% 4. Define amplitude contraint on the SLM (gaussian beam)


% Clear workspace and close all figures
clear all; 
close all;
tic; % Start timing the script

% Add path to the library
addpath('.\lib');

% Propagation parameters
lambda = 1.035e-6; % Wavelength (meters)

%% 1. Define dimensions of the SLM

Nx = 1920; % Number of pixels along x-axis
Ny = 1152; % Number of pixels along y-axis

dx = 9.2e-6; % Pixel size (meters)
dy = 9.2e-6; 
L1 = Nx * dx; % Side size
L2 = Ny * dy;

x = -L1/2 : dx : L1/2 - dx; % Centered screen coordinates in meters
y = -L2/2 : dy : L2/2 - dy; 
[X, Y] = meshgrid(x, y); % Create 2D grid
%|----------- y ^---------------|
%|              |               |
%|--------------|---------------> x
%|              |               |
%|--------------|---------------|


%% 2. Define the region of interest (ROI) for the target image
% here no image is uploaded

Nimg = 60; % Number of pixels in target image. Could be 300
% replaced by [p, q]

% Define centered region of interest (ROI) of target image size
shiftr = 0; % Row shift for ROI
shiftu = 0; % Column shift for ROI
ROIX = Nx/2 - Nimg/2 - shiftu : Nx/2 + Nimg/2 - shiftu; % X-coordinates for ROI
ROIY = Ny/2 - Nimg/2 - shiftr : Ny/2 + Nimg/2 - shiftr; % Y-coordinates for ROI
%ROIX2 = 1920/2 - 1148/2 - 1 : 1920/2 + 1148/2; % previous trials
%ROIY2 = 1152/2 - 1148/2 - 1 : 1152/2 + 1148/2; % previous trials
%^ y[px]------------------------|
%|           ,--------,         |
%|       ROIy|        |         |
%|           '--------'ROIx     |
%|------------------------------> x[px]

% Maks of SLM size, zero outside ROI
mask = zeros(Ny, Nx);
mask(ROIY, ROIX) = 1; 
%^ y[px]------------------------|
%|           ,--------,         |
%| 000...    | 111111 |         |
%|           '--------'  000... |
%|------------------------------> x[px]


%% 3. Amplitude constraint on SLM: gaussian beam
x0 = 0;
y0 = 0;
constant = 1; 
sigma = constant * 2.7 * 10^-3 / 2 * 3; % Standard deviation for Gaussian
Amp=1;

Intensity = (Amp * exp(-(2 * ((X - x0).^2 + (Y - y0).^2 )/ sigma^2))); % Intensity profile
Intensity=Intensity/sum(sum(Intensity));

%% Constrain amplitude and phase in signal region

CircMask(1:Ny, 1:Nx) = 1;
MaskRadius = 5.3e-3;
CircMask(sqrt(X.^2 + Y.^2) >= MaskRadius) = 0; % Circular mask, 1 outside a certain radius

%% 3. Load images or prepare a pattern

%% Create a 10-pixel image
Nimg2 = 10; % Number of pixels in small target image
ROIX3 = Nx/2 - Nimg2/2 - shiftu : Nx/2 + Nimg2/2 - shiftu - 1; % X-coordinates for small image ROI
ROIY3 = Ny/2 - Nimg2/2 - shiftr : Ny/2 + Nimg2/2 - shiftr - 1; % Y-coordinates for small image ROI
mask3 = zeros(Ny, Nx);
mask3(ROIY3, ROIX3) = 1; % Mask for small image ROI

% Load a small target image and scale it
TEN = zeros(Ny, Nx); % Placeholder for 10-pixel image
TEN_p = rgb2gray(imread('pi.png')); % Load 10-pixel target image
TEN(mask3 == 1) = TEN_p; % Apply to mask
TEN = TEN / sum(sum(TEN)); % Normalize intensity
TEN = TEN.^(1/2); % Electric field of 10-pixel image

%% Create helical intensity pattern
% Nxs = 61; Nys = Nxs; % Dimensions of helical pattern
% dxs = 1; % Pixel size
% dys = dxs; % Pixel size
% L1s = Nxs * dxs; % Side length in x
% L2s = Nys * dys; % Side length in y
% xs = -L1s/2 : dxs : L1s/2 - dxs; % Coordinates in x
% ys = -L2s/2 : dys : L2s/2 - dys; % Coordinates in y
% [Xs, Ys] = meshgrid(xs, ys); % Create 2D grid
% helicalpattern = atan2(Ys, Xs); % Generate helical pattern

%% Load and process a sample image (1)
% Target = rgb2gray(imread('Sample_Image.JPG')); % Load target image
% Target = imresize(Target, [Ny Nx]); % Resize to match dimensions
% Target = double(Target); % Convert to double

%% Load and process a sample image (2)
% apple = imread('SinusoidGrating_Period128.bmp'); % apple=rgb2gray(imread('Sample_Image.jpg'));
% apple = imrotate(apple, 270); % Rotate image
% apple = imresize(apple, [Nimg+1 Nimg+1]); % Resize to match target dimensions
% apple = double(apple); % Convert to double
% apple = apple / sum(sum(apple)); % Normalize intensity
% MAX = max(apple(:)); % Get maximum intensity

%% To use when dummy deflection
% maskconstraint_dummy(1:Ny, 1:Nx) = 1;
% maskconstraint_dummy(sqrt((X+(shiftu)*dx).^2+(Y+(shiftr)*dx).^2)>=Nimg/2*dx+0.2*10^-5)=0;
% figure(321);imagesc(maskconstraint_dummy)
% Target(mask==1)=EfieldTarget(Ny/2-Nimg/2:Ny/2+Nimg/2,Nx/2-Nimg/2:Nx/2+Nimg/2);%dummy gaussian
% Target=Target.*maskconstraint_dummy;
% C3=0.45*10^8;%for 300 Nimg
% % C3=2.8*10^8; %for 60 Nimg  % 2.8 ...2mmhole optimized
% C4=(Nimg/2*dx)*1.20;
% iq_phase=C3*(C4-((X+(shiftu)*dx).^2 + (Y+(shiftr)*dx).^2).^(1/2)).^2;

%% Quadratic phase for pondemorotive deflection
% Define mask constraints for ponderomotive deflection
% maskconstraint_gaussian(1:Ny, 1:Nx) = 1;
% maskconstraint_gaussian(sqrt((X + shiftu*dx).^2 + (Y + shiftr*dx).^2) >= Nimg/10*dx + 0.2*10^-5) = 0;

% Target(mask==1)=EfieldTarget(Ny/2-Nimg/2:Ny/2+Nimg/2,Nx/2-Nimg/2:Nx/2+Nimg/2);%gaussian
% Target=Target.*maskconstraint_gaussian;
% C1=220*10^7;
% C2=C1;
% quadraticphase=(C1*(X+(shiftu)*dx).^2 + C2*(Y+(shiftr)*dx).^2);

%% Inverse quadratic phase for ponderomotive deflection
% settings 2: Efield of target (?)
% constant2 = 0.00800; % 0.3 Adjust for different image sizes. Pondemorotive deflection Nimg=60, dummy gaussian Nimg=600
% sigma2 = constant2 * 2.7 * 10^-3; % Standard deviation for target Gaussian
% Amp2 = 1; % Amplitude
% Intensity2 = (Amp2 * exp(-(2 * ((X - x0).^2 + (Y - y0).^2 )/ sigma2^2))); % Intensity profile
% Intensity2 = Intensity2 / sum(sum(Intensity2)); % Normalize intensity
% EfieldTarget = Intensity2.^(1/2); % Electric field
% figure(321321);imagesc(EfieldTarget)

% Define mask constraints for ponderomotive deflection
% maskconstraint_gaussian(1:Ny, 1:Nx) = 1;
% maskconstraint_gaussian(sqrt((X + shiftu*dx).^2 + (Y + shiftr*dx).^2) >= Nimg/10*dx + 0.2*10^-5) = 0;

% Target(mask == 1) = EfieldTarget(Ny/2 - Nimg/2 : Ny/2 + Nimg/2, Nx/2 - Nimg/2 : Nx/2 + Nimg/2); % Gaussian
% Target = Target .* maskconstraint_gaussian; % Apply Gaussian constraints
% C3 = 1 * 10^8; % Quadratic phase coefficient
% C4 = (Nimg/2 * dx) * 0.85; %(Nimg/2*dx)*10; % When no beam block is implemented
% quadraticphase = C3 * (C4 - sqrt((X + shiftu*dx).^2 + (Y + shiftr*dx).^2)).^2; % Compute quadratic phase

%% Quadratic phase for wavefront shaping
% Define mask constraints for wavefront shaping
% maskconstraint(1:Ny, 1:Nx) = 1;
% maskconstraint(sqrt((X + shiftu*dx).^2 + (Y + shiftr*dx).^2) >= Nimg/2*dx + 0.2*10^-5) = 0;

% Target(mask==1)=apple.^(1/2);
% Target=Target.*maskconstraint;
% % C1=4.6*10^7;%for 300 Nimg
% C1=3.5*5*10^7;%for 60 Nimg
% C2=C1;
% quadraticphase=(C1*(X+(shiftu)*dx).^2 + C2*(Y+(shiftr)*dx).^2);

%% Inverse quadratic phase for wavefront shaping (recommended when using beamblock) USE maskconstraint in iterations
% Define mask constraints for wavefront shaping
% maskconstraint(1:Ny, 1:Nx) = 1;
% maskconstraint(sqrt((X + shiftu*dx).^2 + (Y + shiftr*dx).^2) >= Nimg/2*dx + 0.2*10^-5) = 0;

% Target(mask==1)=apple.^(1/2);
% Target=Target.*maskconstraint;
% % C3=0.5*10^8;%for 300 Nimg
% C3=2.8*10^8; %for 60 Nimg  % 2.8 ...2mmhole optimized
% C4=(Nimg/2*dx)*1.20;
% quadraticphase=C3*(C4-((X+(shiftu)*dx).^2 + (Y+(shiftr)*dx).^2).^(1/2)).^2;

%% 4. Initialize phase and target for 10-pixel image

% Electric field constraint
Efield=Intensity.^(1/2);
Efield = Efield .* CircMask; % Apply mask to electric field

% Chose target and its initial phase
Target = TEN; % Set target to 10-pixel image
D = Target.* exp(1i * zeros(Ny, Nx)); % Define initial complex field
% D = Target.*exp(1i*quadraticphase); % waferont shaping
% const=0; % deflection rate
% D = Target.*exp(1i*const); % deflection
% D = Target.*exp(1i*iq_phase);   %dummy deflection

% figure(32);imagesc(angle(D));

% Iteration parameters
k = 1; % Iteration counter, used for quadratic phase initialization and similar
K = 1; % Scaling factor, a support matrix for propagation

error = []; % Initialize error array
iteration_num = 5; % Number of iterations
z = 0.2; % Propagation distance

%% 5. Iterative process for hologram generation
for i = 1: iteration_num
    tic
    % Back-propagation
    A = propagateback(D, -1, lambda, dx, dy, Nx, Ny, X, Y, z, z, z); % Back-propagate
    A = A / sqrt(sum(sum(abs(A).^2))); % Normalize
    as = subplot(3, 5, 1); imagesc(abs(A).^2); title('SLM ideal Int.'); % Display SLM intensity
    as = [as subplot(3, 5, 6)]; imagesc(angle(A)); title('SLM ideal Phase'); % Display SLM phase
    
    % Apply amplitude constraint to the hologram plane
    B = Efield .* exp(1i * angle(A)); % Apply amplitude constraint
    %sum(sum(abs(B).^2))
    as = [as subplot(3, 5, 2)]; imagesc(abs(B).^2); title('SLM true Int.'); % Display true intensity
    as = [as subplot(3, 5, 7)]; imagesc(angle(B)); title('SLM true Phase'); % Display true phase
    
    % Forward-propagation
    %C = 1/sqrt(Nx*Ny)*fftshift(fft2(fftshift(B)));
    C = propagateforward(B, 1, lambda, dx, dy, Nx, Ny, X, Y, z, z, z); % Forward-propagate
    C = C / sqrt(sum(sum(abs(C).^2))); % Normalize
    as = [as subplot(3, 5, 3)]; imagesc(abs(C).^2); title('Image Int.'); % Display image intensity
    as = [as subplot(3, 5, 8)]; imagesc(angle(C)); title('Image Phase'); % Display image phase
    
    D = C; % Copy current field for next iteration
    
    % Apply constraints based on image type
        % double constraint
    %D(maskconstraint_gaussian==1) = abs(Target(maskconstraint_gaussian==1)).*exp(1i*k*quadraticphase(maskconstraint_gaussian==1));%signal window
    %D(mask3==1) = abs(Target(mask3==1)).*exp(1i*k*Phase_TEN(mask3==1));%signal window TEN pixels img
        % single constraint
    %D(maskconstraint==1) = K*abs(Target(maskconstraint==1)).*exp(1i*angle(C(maskconstraint==1)));%wavefront shaping
    %D(maskconstraint_gaussian==1) = K*abs(Target(maskconstraint_gaussian==1)).*exp(1i*angle(C(maskconstraint_gaussian==1)));%deflection
    D(mask3 == 1) = K * abs(Target(mask3 == 1)) .* exp(1i * angle(C(mask3 == 1))); % 10-pixel image constraint
    %D(maskcons1traint_dummy==1) = K*abs(Target(maskconstraint_dummy==1)).*exp(1i*angle(C(maskconstraint_dummy==1)));%dummy
    %D = K*abs(Target).*exp(1i*angle(C));
    %D=D/sqrt(sum(sum(abs(D).^2)));

    as = [as subplot(3, 5, 4)]; imagesc(abs(D).^2); title('Image ideal Int.'); % Display ideal intensity
    as = [as subplot(3, 5, 9)]; imagesc(angle(D)); title('Image enforced phase'); % Display enforced phase
    
    % Calculate error and efficiency
    error = sum(sum(abs(abs(C .* mask).^2 - abs(Target .* mask).^2))); % Compute error
    %correlationCODE=corr2(abs(C.*mask).^2,abs(Target.*mask).^2);
    deff = sum(sum(abs(C .* mask).^2)) / sum(sum(abs(C).^2)); % Compute efficiency
    subplot(3, 5, 5); hold on; plot(i, error, 'bo'); hold off; title('Error'); % Plot error
    subplot(3, 5, 10); hold on; plot(i, deff, 'bo'); hold off; title('Efficiency'); % Plot efficiency
    %subplot(3,5,11); hold on ; plot(i,correlation,'bo'); hold off; title('Correlation');

    toc % End timing for this iteration
    drawnow; % Update figures
    linkaxes(as); % Link axes for simultaneous zoom
end

% Display final image
figure(3211); imagesc(abs(C).^2)
