function yy = ylm1(L, m , theta0, phi0) 
%% fully normalized spherical harmonics 
% y = ylm(L, m , theta, phi) 
%                   [ (2*L+1) (L-m)!   ]
% ylm = SQRT |---------------|  Plm(cos\theta) exp(I*m*phi) 
%                   [   4 *pi  (L+m)!   ]
% this definition is consistent with Mathemaitca(R) definition 
% 
% yingcai zheng, 2014 
% yc.zheng@gmail.com 
% (copyrighted) 
% to calculate the ylm of a series of theta phi
% modified by Yuan Tian
nth = numel(theta0); 
nph = numel(phi0); 

theta = reshape(theta0, 1, nth ); 
phi = reshape(phi0, 1, nph); 
% theta=theta0;
% phi=phi0;
mabs = abs(m) ; 
if mabs > L 
    disp ( 'wrong ') ; 
    pause(100000); 
    return
end
y = legendre(L, cos(theta), 'norm'); 
y = (-1)^mabs* y(mabs+1,:) /sqrt(2*pi); % get m-th order
if m < 0
    y = (-1)^mabs *y; 
end

im = complex(0,1);
y2 = exp(im*m*phi) ;
yy = (y).* y2; 
% yy=transpose(yy);
end
