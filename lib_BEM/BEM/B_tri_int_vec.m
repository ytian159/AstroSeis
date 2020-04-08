function [ST]=B_tri_int_vec(vp,w,...
    xi,yi,zi,wi,xs,ys,zs,area,STs,k)

Gpn = greens_func_p (vp,w,...
    xi,yi,zi,xs,ys,zs);
fac = 2; 
ST = sum(Gpn.*wi,2) .*area*fac;

ST(k)=STs;

end