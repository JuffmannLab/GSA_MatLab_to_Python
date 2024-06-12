clear all; close all;
tic;

lambda =1.035e-6;

%% specify SLM dimensions
Nx=1920;
Ny=1152;
dx=(9.2e-6);
dy=(9.2e-6);
L1=1920*dx; %side length
L2=1152*dy;
x=-L1/2:dx:L1/2-dx; %src coords
y=-L2/2:dy:L2/2-dy;
[X,Y] = meshgrid(x,y);

%% Load and prepare target image
%Target = rgb2gray(imread('Smiley50px.png')); 
Target = rgb2gray(imread('2gauss.png')); 
[p,q]=size(Target);
shiftr=0;
shiftu=0;
ROIX=Nx/2-q/2-shiftu:Nx/2+q/2-shiftu-1;
ROIY=Ny/2-p/2-shiftr:Ny/2+p/2-shiftr-1;
Targetn=zeros(Ny,Nx);
Targetn(Ny/2-p/2:Ny/2+p/2-1,Nx/2-q/2:Nx/2+q/2-1)=abs(255-Target);
%Targetn(ROIX,ROIY)=Target;
Targetn = Targetn/sum(sum(Targetn));
Targetn = Targetn.^(1/2);

%% Define Laser field
x0 = 0;
y0 = 0;
sigma1=2.7*10^-3/2*3;
sigma2=2.7*10^-3/2*3;

Amp=1;
res = (2*(X-x0).^2./(sigma1^2) + 2*(Y-y0).^2./(sigma1^2));
Intensity = ((Amp  * exp(-res)));
Intensity=Intensity/sum(sum(Intensity));
Efield=Intensity.^(1/2);

mask=zeros(Ny,Nx);
mask(ROIY,ROIX)=1;

%constrain the amplitude and phase only in the signal region
SLM=ones(1152,1920);
Tf(1:Ny,1:Nx)=1;
Tf(sqrt(X.^2+Y.^2)>=5.3e-3)=0;

Phase_TEN=zeros(Ny,Nx);
D = Targetn.*exp(1i*Phase_TEN); 
Efield=Efield.*Tf;

error = [];
iteration_num= 20;

z = 0.1;

for i=1:iteration_num
tic
       
        A = propagateback(D,-1,lambda,dx,dy,Nx,Ny,X,Y,z,z,z);
        A=A./sqrt(sum(sum(abs(A).^2)));
        as = subplot(3,5,1); imagesc(abs(A).^2); title('SLM ideal Int.'); 
        as = [as subplot(3,5,6)]; imagesc(angle(A)); title('SLM ideal Phase'); 

        B = Efield.*exp(1i*angle(A));%apply amplitude constraint to the hologram plane
        %sum(sum(abs(B).^2))
        as = [as subplot(3,5,2)]; imagesc(abs(B).^2); title('SLM true Int.'); 
        as = [as subplot(3,5,7)]; imagesc(angle(B)); title('SLM true Phase'); 
% 
        %C = 1/sqrt(Nx*Ny)*fftshift(fft2(fftshift(B)));
        C = propagateforward(B,1,lambda,dx,dy,Nx,Ny,X,Y,z,z,z);
        C = C./sqrt(sum(sum(abs(C).^2)));
        as = [as subplot(3,5,3)]; imagesc(abs(C).^2); title('Image Int.'); 
        as = [as subplot(3,5,8)]; imagesc(angle(C)); title('Image Phase'); 

        D = Targetn.*exp(1i*angle(C)); %noise window
      
        as = [as subplot(3,5,4)]; imagesc(abs(D).^2); title('Image ideal Int.');
        as = [as subplot(3,5,9)]; imagesc(angle(D)); title('Image enforced phase');

        error = sum(sum(abs(abs(C.*mask).^2 - abs(Targetn.*mask).^2)));

        deff = sum(sum(abs(C.*mask).^2))/sum(sum(abs(C).^2));
        subplot(3,5,5); hold on ; plot(i,error,'bo'); hold off; title('Error');
        subplot(3,5,10); hold on ; plot(i,deff,'bo'); hold off; title('Efficiency');

        toc
        drawnow;
        linkaxes(as);
     
    end

    figure(3211); imagesc(abs(C).^2)
    figure(1111); imagesc(angle(C))

%% Image to display on the slm
data=angle(C);
% Normalize the data to the range [0, 1]
min_val = min(data(:));
max_val = max(data(:));
normalized_data = (data - min_val) / (max_val - min_val);

% Convert normalized data to 8-bit unsigned integers
image_data = uint8(normalized_data * 255);

% Save the image as BMP
imwrite(image_data, '2gauss-marius.bmp');


%%
function[u2]=propagate(u1,sgn,lambda,dx,dy,Nx,Ny,X,Y,z1,z2,f)
% z1=20e-2; % Distance SLM to Lens
% z2=20e-2; % Distance Lens to screen
% f=20e-2;  % focal length of lens
SLM=ones(1152,1920);
[nor,noc]=size(SLM);
Tf(1:nor,1:noc)=1;
Tf(sqrt(X.^2+Y.^2)>=10e-3)=0; %pupil function of lens
% Th(1:nor,1:noc)=1;
% Th(sqrt(X.^2+Y.^2)<=8e-4)=0; %pupil function of mirror
k=2*pi/lambda;

u_dummy=propTF2(u1,dx*Nx,dy*Ny,lambda,sgn*z1);
u_dummy=u_dummy.*Tf.*exp(-sgn*1i*k/2/f*(X.^2+Y.^2));
 u2=propTF2(u_dummy,dx*Nx,dy*Ny,lambda,sgn*z2);
%  u2=propTF1(u_dummy,dx*Nx,dy*Ny,lambda,sgn*z2/2);
%  u2=u2.*Th;
%  u2=propTF1(u2,dx*Nx,dy*Ny,lambda,sgn*z2/2);
end
function[u2]=propagateback(u1,sgn,lambda,dx,dy,Nx,Ny,X,Y,z1,z2,f)
% z1=20e-2; % Distance SLM to Lens
% z2=20e-2; % Distance Lens to screen
% f=20e-2;  % focal length of lens
shiftr=-0;
shiftu=0;
SLM=ones(1152,1920);
[nor,noc]=size(SLM);
Tf(1:nor,1:noc)=1;
Tf(sqrt(X.^2+Y.^2)>=12.5e-3)=0; %pupil function of lens
Th(1:nor,1:noc)=1;
% Th(sqrt((X+(shiftu)*dx).^2+(Y+(shiftr)*dx).^2)<=1e-3)=0; 
k=2*pi/lambda;


 u2=propTF2(u1,dx*Nx,dy*Ny,lambda,sgn*z2/2);
 u2=u2.*Th;
 u2=propTF2(u2,dx*Nx,dy*Ny,lambda,sgn*z2/2);

 u2=u2.*Tf.*exp(-sgn*1i*k/2/f*(X.^2+Y.^2));
 u2=propTF2(u2,dx*Nx,dy*Ny,lambda,sgn*z2);

end
function[u2]=propagateforward(u1,sgn,lambda,dx,dy,Nx,Ny,X,Y,z1,z2,f)
% z1=20e-2; % Distance SLM to Lens
% z2=20e-2; % Distance Lens to screen
% f=20e-2;  % focal length of lens
shiftr=-0;
shiftu=0;
SLM=ones(1152,1920);
[nor,noc]=size(SLM);
Tf(1:nor,1:noc)=1;
Tf(sqrt(X.^2+Y.^2)>=12.5e-3)=0; %pupil function of lens
Th(1:nor,1:noc)=1;
% Th(sqrt((X+(shiftu)*dx).^2+(Y+(shiftr)*dx).^2)<=1e-3)=0; 
k=2*pi/lambda;

u_dummy=propTF2(u1,dx*Nx,dy*Ny,lambda,sgn*z1);
u_dummy=u_dummy.*Tf.*exp(-sgn*1i*k/2/f*(X.^2+Y.^2));

 u2=propTF2(u_dummy,dx*Nx,dy*Ny,lambda,sgn*z2/2);
 u2=u2.*Th;
 u2=propTF2(u2,dx*Nx,dy*Ny,lambda,sgn*z2/2);
end

function[u2]=propagatetohole(u1,sgn,lambda,dx,dy,Nx,Ny,X,Y,z1,z2,f)
% z1=20e-2; % Distance SLM to Lens
% z2=20e-2; % Distance Lens to screen
% f=20e-2;  % focal length of lens
shiftr=150;
SLM=ones(1152,1920);
[nor,noc]=size(SLM);
Tf(1:nor,1:noc)=1;
Tf(sqrt(X.^2+Y.^2)>=12.5e-3)=0; %pupil function of lens
Th(1:nor,1:noc)=1;
Th(sqrt(X.^2+(Y+(shiftr)*dx).^2)<=0.5e-3)=0; 
k=2*pi/lambda;

u_dummy=propTF2(u1,dx*Nx,dy*Ny,lambda,sgn*z1);
u_dummy=u_dummy.*Tf.*exp(-sgn*1i*k/2/f*(X.^2+Y.^2));

 u2=propTF2(u_dummy,dx*Nx,dy*Ny,lambda,sgn*z2/2);
%  u2=u2.*Th;
%  u2=propTF1(u2,dx*Nx,dy*Ny,lambda,sgn*z2/2);
end