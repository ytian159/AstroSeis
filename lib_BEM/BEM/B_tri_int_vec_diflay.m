function [ST]=B_tri_int_vec_diflay(vp,w,...
    xi,yi,zi,wi,xs,ys,zs,area)

Gpn = greens_func_p (vp,w,...
    xi,yi,zi,xs,ys,zs);
fac = 2; 
ST = sum(Gpn.*wi,2) .*area*fac;



end