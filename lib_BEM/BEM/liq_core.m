function [A,b]=liq_core(face1,face2,w,lamda1,mu1,rho1,lamda2,mu2,rho2,Q,...
    nint,nxi,P0,u01,u02,Smat,w0)
face2w1=face2;
for i=1:length(face2)
    face2w1(i).nvec=-face2(i).nvec;
end
%% interaction matrix
Tst_st11=cal_T_st(face1,face1,w,lamda1,mu1,rho1,Q,nint,nxi);
Tst_st22=cal_T_st(face2w1,face2w1,w,lamda1,mu1,rho1,Q,nint,nxi);
T_st12=cal_T_st(face1,face2w1,w,lamda1,mu1,rho1,Q,nint,nxi);
T_st21=cal_T_st(face2w1,face1,w,lamda1,mu1,rho1,Q,nint,nxi);

%G_st11=cal_G_st(face1,face1,w,lamda,mu,rho,Q,nint,nxi);
G_st12=cal_G_st(face1,face2w1,w,lamda1,mu1,rho1,Q,nint,nxi);
G_st22=cal_G_st(face2w1,face2w1,w,lamda1,mu1,rho1,Q,nint,nxi);

Ast_st22=cal_A_st(face2,face2,w,lamda2,mu2,rho2,Q,nint,nxi);

B_st22=cal_B_st(face2,face2,w,lamda2,mu2,rho2,Q,nint,nxi);
%% 
sclae_fac=(rho2*sqrt(lamda1+2*mu1)/sqrt(rho2)*w0);
zrmat=zeros(length(face2),length(face1)*3);
%%
A=[Ast_st22/sclae_fac,-rho2*w^2*B_st22*Smat/sclae_fac,zrmat;
    -G_st22*Smat.',Tst_st22,T_st21;
    -G_st12*Smat.',T_st12,Tst_st11];
b=[P0/sclae_fac;u02;u01];
end