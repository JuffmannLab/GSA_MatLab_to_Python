clear all; close all;
tic;

% Add path to my library
addpath('C:\Users\thomas\Desktop\temp_gsa\lib')

% test=imread('laserimage_outcoupled.bmp');

Nx=1920;
Ny=1152;

% target image roi
% Nimg=300; %nr of pixels in target image
Nimg=60; %nr of pixels in target image
shiftr=0;
shiftu=0;
ROIX=Nx/2-Nimg/2-shiftu:Nx/2+Nimg/2-shiftu;
ROIY=Ny/2-Nimg/2-shiftr:Ny/2+Nimg/2-shiftr;


ROIX2=1920/2-1148/2-1:1920/2+1148/2;
ROIY2=1152/2-1148/2-1:1152/2+1148/2;


%%


mask=zeros(Ny,Nx);
mask(ROIY,ROIX)=1;

mask2=zeros(Ny,Nx);
mask2(ROIY2,ROIX2)=1;

lambda =1.035e-6;

dx=(9.2e-6);
dy=(9.2e-6);
L1=1920*dx; %side length
L2=1152*dy;

x=-L1/2:dx:L1/2-dx; %src coords
y=-L2/2:dy:L2/2-dy;
[X,Y] = meshgrid(x,y);
x0 = 0;
y0 = 0;

constant=1;
 sigma1=constant*2.7*10^-3/2*3;
 sigma2=constant*2.7*10^-3/2*3;

Amp=1;
res = (2*(X-x0).^2./(sigma1^2) + 2*(Y-y0).^2./(sigma1^2));
Intensity = ((Amp  * exp(-res)));
Intensity=Intensity/sum(sum(Intensity));
Efield=Intensity.^(1/2);



%% define target gaussian
 constant=0.00800; %pondemorotive deflection Nimg=60
%  constant=0.3; %dummy gaussian Nimg=600
 sigma1=constant*2.7*10^-3;
 res = (2*(X-x0).^2./(sigma1^2) + 2*(Y-y0).^2./(sigma1^2));
 Intensity = ((Amp  * exp(-res)));
 Intensity=Intensity/sum(sum(Intensity));
 EfieldTarget=Intensity.^(1/2);
% figure(321321);imagesc(EfieldTarget)


%%

% define noise region a.k.a zeros
Target=rgb2gray(imread('Sample_Image.JPG'));
Target=imresize(Target, [Ny Nx]);
Target=double(Target);
Target=zeros*(Target);
Zernikepadded=Target;
Zernikepadded_tilt=Target;
Zernikepadded_defocus=Target;
Zernikepadded_spherical=Target;
TEN = Target;

%% create 10 pixels images
Nimg2=10; %nr of pixels in target image
shiftr=0;
shiftu=0;
ROIX3=Nx/2-Nimg2/2-shiftu:Nx/2+Nimg2/2-shiftu-1;
ROIY3=Ny/2-Nimg2/2-shiftr:Ny/2+Nimg2/2-shiftr-1;
mask3=zeros(Ny,Nx);
mask3(ROIY3,ROIX3)=1;


TEN_p = rgb2gray(imread('pi.png')); %Targetimage with 10 pixels
TEN(mask3==1) = TEN_p;
TEN = TEN/sum(sum(TEN));
TEN = TEN.^(1/2);
%% 
%% helical intensity pattern
Nxs=61; Nys=Nxs;
dxs=1;
dys=dxs;
L1s=Nxs*dxs; %side length
L2s=Nys*dys;
xs=-L1s/2:dxs:L1s/2-dxs;
ys=-L2s/2:dys:L2s/2-dys;
[Xs,Ys]=meshgrid(xs,ys);
helicalpattern=atan2(Ys,Xs);

%%
% apple=rgb2gray(imread('Sample_Image.jpg'));
apple=imread('SinusoidGrating_Period128.bmp');
% apple=imread('z.JPG');
apple= imrotate(apple,270);
apple=imresize(apple, [Nimg+1 Nimg+1]);
apple=double(apple);
apple=apple/sum(sum(apple));
MAX=max(apple(:));


%constrain the amplitude and phase only in the signal region

SLM=ones(1152,1920);
[nor,noc]=size(SLM);
Tf2(1:nor,1:noc)=0;
Tf2(sqrt(X.^2+Y.^2)>=5.3e-3)=1;

Tf(1:nor,1:noc)=1;
maskconstraint=Tf;
maskconstraint_gaussian=Tf;
maskconstraint_dummy=Tf;
Tf(sqrt(X.^2+Y.^2)>=5.3e-3)=0;

maskconstraint(sqrt((X+(shiftu)*dx).^2+(Y+(shiftr)*dx).^2)>=Nimg/2*dx+0.2*10^-5)=0; %for wavefront shaping


%% to use when ponderomotive deflection
maskconstraint_gaussian(sqrt((X+(shiftu)*dx).^2+(Y+(shiftr)*dx).^2)>=Nimg/10*dx+0.2*10^-5)=0;

%% to use when dummy deflection
% maskconstraint_dummy(sqrt((X+(shiftu)*dx).^2+(Y+(shiftr)*dx).^2)>=Nimg/2*dx+0.2*10^-5)=0;
% figure(321);imagesc(maskconstraint_dummy)
% Target(mask==1)=EfieldTarget(Ny/2-Nimg/2:Ny/2+Nimg/2,Nx/2-Nimg/2:Nx/2+Nimg/2);%dummy gaussian
% Target=Target.*maskconstraint_dummy;
% C3=0.45*10^8;%for 300 Nimg
% % C3=2.8*10^8; %for 60 Nimg  % 2.8 ...2mmhole optimized
% C4=(Nimg/2*dx)*1.20;
% iq_phase=C3*(C4-((X+(shiftu)*dx).^2 + (Y+(shiftr)*dx).^2).^(1/2)).^2;

%% quadratic phase for pondemorotive deflection
% Target(mask==1)=EfieldTarget(Ny/2-Nimg/2:Ny/2+Nimg/2,Nx/2-Nimg/2:Nx/2+Nimg/2);%gaussian
% Target=Target.*maskconstraint_gaussian;
% C1=220*10^7;
% C2=C1;
% quadraticphase=(C1*(X+(shiftu)*dx).^2 + C2*(Y+(shiftr)*dx).^2);

%% inverse quadratic phase for pondemorotive deflection (recommended when using beamblock) USE maskconstraint_gaussian in iterations
Target(mask==1)=EfieldTarget(Ny/2-Nimg/2:Ny/2+Nimg/2,Nx/2-Nimg/2:Nx/2+Nimg/2);%gaussian
Target=Target.*maskconstraint_gaussian;
C3=1*10^8;
C4=(Nimg/2*dx)*0.85; %when no bb implemented
% C4=(Nimg/2*dx)*10; %when no bb implemented
quadraticphase=C3*(C4-((X+(shiftu)*dx).^2 + (Y+(shiftr)*dx).^2).^(1/2)).^2;

%% quadratic phase for wavefront shaping
% Target(mask==1)=apple.^(1/2);
% Target=Target.*maskconstraint;
% % C1=4.6*10^7;%for 300 Nimg
% C1=3.5*5*10^7;%for 60 Nimg
% C2=C1;
% quadraticphase=(C1*(X+(shiftu)*dx).^2 + C2*(Y+(shiftr)*dx).^2);


%% inverse quadratic phase for wavefront shaping (recommended when using beamblock) USE maskconstraint in iterations
% Target(mask==1)=apple.^(1/2);
% Target=Target.*maskconstraint;
% % C3=0.5*10^8;%for 300 Nimg
% C3=2.8*10^8; %for 60 Nimg  % 2.8 ...2mmhole optimized
% C4=(Nimg/2*dx)*1.20;
% quadraticphase=C3*(C4-((X+(shiftu)*dx).^2 + (Y+(shiftr)*dx).^2).^(1/2)).^2;

%%
const=0;
% D = Target.*exp(1i*quadraticphase); %waferont shaping
% D = Target.*exp(1i*const); %deflection
%% %10 pixels image
Phase_TEN=zeros(Ny,Nx);
Target = TEN;
D = TEN.*exp(1i*Phase_TEN); 
%%
% D = Target.*exp(1i*iq_phase);   %dummy deflection
% figure(32);imagesc(angle(D));
Efield=Efield.*Tf;


k=1;
K=1;

error = [];
iteration_num= 10;
clims = [0 MAX];

for i=1:iteration_num
tic
        %A = sqrt(Nx*Ny)*ifftshift(ifft2(ifftshift(D)));
        A = propagateback(D,-1,lambda,dx,dy,Nx,Ny,X,Y,0.2,0.2,0.2);
        A=A./sqrt(sum(sum(abs(A).^2)));
        as = subplot(3,5,1); imagesc(abs(A).^2); title('SLM ideal Int.'); 
        as = [as subplot(3,5,6)]; imagesc(angle(A)); title('SLM ideal Phase'); 

        B = Efield.*exp(1i*angle(A));%apply amplitude constraint to the hologram plane
        %sum(sum(abs(B).^2))
        as = [as subplot(3,5,2)]; imagesc(abs(B).^2); title('SLM true Int.'); 
        as = [as subplot(3,5,7)]; imagesc(angle(B)); title('SLM true Phase'); 
% 
        %C = 1/sqrt(Nx*Ny)*fftshift(fft2(fftshift(B)));
        C = propagateforward(B,1,lambda,dx,dy,Nx,Ny,X,Y,0.2,0.2,0.2);
        C = C./sqrt(sum(sum(abs(C).^2)));
        as = [as subplot(3,5,3)]; imagesc(abs(C).^2); title('Image Int.'); 
        as = [as subplot(3,5,8)]; imagesc(angle(C)); title('Image Phase'); 

        D = C; %noise window
        %% double constraint
%          D(maskconstraint_gaussian==1) = abs(Target(maskconstraint_gaussian==1)).*exp(1i*k*quadraticphase(maskconstraint_gaussian==1));%signal window
%            D(mask3==1) = abs(Target(mask3==1)).*exp(1i*k*Phase_TEN(mask3==1));%signal window TEN pixels img
        %% single constraint
%         D(maskconstraint==1) = K*abs(Target(maskconstraint==1)).*exp(1i*angle(C(maskconstraint==1)));%wavefront shaping
%         D(maskconstraint_gaussian==1) = K*abs(Target(maskconstraint_gaussian==1)).*exp(1i*angle(C(maskconstraint_gaussian==1)));%deflection
        D(mask3==1) = K*abs(Target(mask3==1)).*exp(1i*angle(C(mask3==1)));%TEN pixels img
%         D(maskcons1traint_dummy==1) = K*abs(Target(maskconstraint_dummy==1)).*exp(1i*angle(C(maskconstraint_dummy==1)));%dummy
%         D = K*abs(Target).*exp(1i*angle(C)); %GSA

        %D=D/sqrt(sum(sum(abs(D).^2)));
        as = [as subplot(3,5,4)]; imagesc(abs(D).^2); title('Image ideal Int.');
        as = [as subplot(3,5,9)]; imagesc(angle(D)); title('Image enforced phase');

        error = sum(sum(abs(abs(C.*mask).^2 - abs(Target.*mask).^2)));
%         correlationCODE=corr2(abs(C.*mask).^2,abs(Target.*mask).^2);
        deff = sum(sum(abs(C.*mask).^2))/sum(sum(abs(C).^2));
        subplot(3,5,5); hold on ; plot(i,error,'bo'); hold off; title('Error');
        subplot(3,5,10); hold on ; plot(i,deff,'bo'); hold off; title('Efficiency');
%         subplot(3,5,11); hold on ; plot(i,correlation,'bo'); hold off; title('Correlation');

        toc
        drawnow;
        linkaxes(as);
      
    end

    figure(3211); imagesc(abs(C).^2)

