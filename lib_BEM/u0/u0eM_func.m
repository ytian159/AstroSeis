function [u0]=u0eM_func(face,w,rho,mu,lamda,xs,ys,zs,Q,M)
numface=length(face);
u0=zeros(numface*3,1);
vp0=sqrt((lamda+2*mu)/rho);
vs0=sqrt(mu/rho);
Qp=2.5*Q;
vp=vp0/(1+complex(0,1)*0.5/Qp);
vs=vs0/(1+complex(0,1)*0.5/Q);
for j =1: numface
    x1 = face(j).ic(1); x2 = face(j).ic(2); x3=face(j).ic(3);
    [gi1,gi2,gi3]=greens_deri_src(vp,vs,rho,w,...
    x1,x2,x3,xs,ys,zs);
    u0(j)=sum(sum(gi1.*M));
    u0(j+numface)=sum(sum(gi2.*M));
    u0(j+numface*2)=sum(sum(gi3.*M));
end
%%