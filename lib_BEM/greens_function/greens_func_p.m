function [Ep]=greens_func_p(vp,w,x1,x2,x3,xs,ys,zs) 
im = complex(0,1); 
r = sqrt( (x1-xs).^2 + (x2-ys).^2 + (x3-zs).^2); 
kp = w./vp; 
Ep = exp(im.*kp.*r)./(4.*pi.*r);
end