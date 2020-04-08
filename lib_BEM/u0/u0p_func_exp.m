function [u0]=u0p_func_exp(face,w,vp0,xs,ys,zs,Q)
numface=length(face);
Qp=2.5*Q;
vp=vp0/(1+complex(0,1)*0.5/Qp);

for j =1: numface
    x1 = face(j).ic(1); x2 = face(j).ic(2); x3=face(j).ic(3);
    
    %[ux,uy,uz] = disp_due_to_explosion (vp,vs,Q, rho,w,x1,x2,x3,xsource,ysource,zsource);
    ur0=greens_func_p(vp,w,x1,x2,x3,xs,ys,zs);
    %u0(3*j-2:3*j) =[ux uy uz];
    u0(j)=ur0;   
end
u0=u0(:);

end