function [u0]=u0e_func(face,w,rho,mu,lamda,xs,ys,zs,Q,fsrc)
numface=length(face);
u0=zeros(numface*3,1);
for j =1: numface
    x1 = face(j).ic(1); x2 = face(j).ic(2); x3=face(j).ic(3);
    %[ux,uy,uz] = disp_due_to_explosion (vp,vs,Q, rho,w,x1,x2,x3,xsource,ysource,zsource);
    Grs=greens_functionQ(w,rho,mu,lamda,[xs,ys,zs],x1,x2,x3,Q);
    ur0=Grs*fsrc.';
    %u0(3*j-2:3*j) =[ux uy uz];
    u0(j)=ur0(1);
    u0(j+numface)=ur0(2);
    u0(j+numface*2)=ur0(3);
    
end
u0=u0(:);

end