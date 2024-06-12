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