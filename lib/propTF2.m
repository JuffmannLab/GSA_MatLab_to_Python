 function[u2]=propTF1(u1,L1,L2,lambda,z);
 % propagation - transfer function approach
 % assumes same x and y side lengths and
 % uniform sampling
 % u1 - source plane field
 % L - source and observation plane side length
 % lambda - wavelength
 % z - propagation distance
 % u2 - observation plane field
 %  dx=L/M; sample interval

 %  [M,N]=size(u1); %get input field array size
 
 dx=(9.2e-6);
 dy=(9.2e-6); %sample interval

 %L1=1152*(9.2e-6); %side length
 %L2=1152*(9.2e-6);
 L1=1920*dx;     %side length
 L2=1152*dx;

 
 k=2*pi/lambda; %wavenumber

 fx=-1/(2*dx):1/L1:1/(2*dx)-1/L1; %freq coords
 fy=-1/(2*dy):1/L2:1/(2*dy)-1/L2;
 [FX,FY]=meshgrid(fx,fy);

 H=exp(-j*pi*lambda*z*(FX.^2+FY.^2)); %trans func
 H=fftshift(H); %shift trans func
 U1=fft2(fftshift(u1)); %shift, fft src field
 U2=H.*U1; %multiply
 u2=ifftshift(ifft2(U2)); %inv fft, center obs field

 end