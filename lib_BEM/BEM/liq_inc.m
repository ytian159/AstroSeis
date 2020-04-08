function [A,b]=liq_inc(face,w,lamda1,mu1,rho1,lamda2,mu2,rho2,Q,nint,nxi,Smat,P0,u0,w0)
facew1=face;

for i=1:length(face)
    facew1(i).nvec=-face(i).nvec; %flip the surface normal direction for CMB
end

Tst_st22=cal_T_st(facew1,facew1,w,lamda1,mu1,rho1,Q,nint,nxi);

Ast_st22=cal_A_st(face,face,w,lamda2,mu2,rho2,Q,nint,nxi);

G_st22=cal_G_st(facew1,facew1,w,lamda1,mu1,rho1,Q,nint,nxi);

B_st22=cal_B_st(face,face,w,lamda2,mu2,rho2,Q,nint,nxi);
%%
sclae_fac=(rho2*sqrt(lamda1+2*mu1)/sqrt(rho2)*w0); % s=rho c w0
A=[Ast_st22/sclae_fac,-rho2*w^2*B_st22*Smat/sclae_fac;-G_st22*Smat.',Tst_st22];
b=[P0/sclae_fac;u0];

end