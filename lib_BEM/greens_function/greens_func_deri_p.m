function [Gpn]=greens_func_deri_p(vp,w,x1,x2,x3,xs,ys,zs,n1,n2,n3) 
im = complex(0,1); 
r = sqrt( (x1-xs).^2 + (x2-ys).^2 + (x3-zs).^2); 
kp = w./vp; 
Ep = exp(im*kp*r)./(4.*pi.*r);
g1 = (x1-xs)./r; 
g2 = (x2-ys)./r; 
g3 = (x3-zs)./r; 
Gpn=Ep.*(-1./r+im*kp).*(g1.*n1+g2.*n2+g3.*n3);
end