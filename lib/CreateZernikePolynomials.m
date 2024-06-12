function Zernike=CreateZernikePolynomials(dx,k)
      
      n = [0  1  1  2  2  2  3  3  3  3  4  4  4  4  4  5  5  5  5  5  5  6  6  6  6  6  6  6];
      m = [0 -1  1 -2  0  2 -3 -1  1  3 -4 -2  0  2  4 -5 -3 -1  1  3  5 -6 -4 -2  0  2  4  6];
      
%       k = 5;
%       dx=0.00174;
      x1=-1:dx:1;
      [X1,Y1] = meshgrid(x1,x1);
      [theta,r] = cart2pol(X1,Y1);
      idx = r<=1;
      z = nan(size(X1));
      y = zernfun(n,m,r(idx),theta(idx));
    
      z(idx)=y(:,k);
      Zernike=z*pi;
%       figure();imagesc(Zernike); axis off;
%       pause(2);
     
%       z(isnan(z))=0;
      
%% 
%       Mat=-(z+pi)./2/pi;
%        imwrite(Mat, 'Mat-.bmp')
      
      %make zernike in grayvalues
%       Zernike = z - min(z(:));
%       Zernike = Zernike ./ max(Zernike(:));
%       Zernike(isnan(Zernike))=0;
%       figure(); imagesc(x1,x1,Zernike);
     
%       imwrite(Zernike,'Zernike.bmp');
end
