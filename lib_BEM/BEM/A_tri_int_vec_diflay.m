function [ST]=A_tri_int_vec_diflay(vp,w,...
    xi,yi,zi,wi,xs,ys,zs,area,n1,n2,n3,nsub)
n1=repmat(n1,1,nsub);
n2=repmat(n2,1,nsub);
n3=repmat(n3,1,nsub);
Gpn = greens_func_deri_p (vp,w,...
    xi,yi,zi,xs,ys,zs,n1,n2,n3);
fac = 2; 
ST = sum(Gpn.*wi,2) .*area.*fac;

%ST(k)=STs;
end