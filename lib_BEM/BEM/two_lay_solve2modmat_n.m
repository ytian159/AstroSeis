function [u1,Tst_st11,P]=two_lay_solve2modmat_n(face1,face2,w,lamda1,mu1,rho1,lamda2,mu2,rho2,Q,nint,nxi,Sst,P0,u01,u02,Smat)
    face2w1=face2;
    for i=1:length(face2)
        face2w1(i).nvec=-face2(i).nvec;
    end
    Tst_st11=cal_T_st(face1,face1,w,lamda1,mu1,rho1,Q,nint,nxi);
    Tst_st22=cal_T_st(face2w1,face2w1,w,lamda1,mu1,rho1,Q,nint,nxi);
    T_st12=cal_T_st(face1,face2w1,w,lamda1,mu1,rho1,Q,nint,nxi);
    T_st21=cal_T_st(face2w1,face1,w,lamda1,mu1,rho1,Q,nint,nxi);
    
    %G_st11=cal_G_st(face1,face1,w,lamda,mu,rho,Q,nint,nxi);
    G_st12=cal_G_st(face1,face2w1,w,lamda1,mu1,rho1,Q,nint,nxi);
    G_st22=cal_G_st(face2w1,face2w1,w,lamda1,mu1,rho1,Q,nint,nxi);
    %%
    Ast_st22=cal_A_st(face2,face2,w,lamda2,mu2,rho2,Q,nint,nxi);
    %%
    B_st22=cal_B_st(face2,face2,w,lamda2,mu2,rho2,Q,nint,nxi);
    %%
    %Sst2=inv(damping_factor*eye(nS)+SS)*Smat.';
    Bst=1/(rho1*w^2)*Sst*B_st22^-1;
    BAst=Bst*Ast_st22;
    %%
    P0st=Bst*P0;
    u01T12P0=u01+T_st12*P0st;
    u02T22P0=u02+Tst_st22*P0st;
%     T12G12=T_st12*BAst+G_st12*Smat.'; % question here: might be T_st12*BAst-G_st12*Smat.' here
%     T22G22=Tst_st22*BAst+G_st22*Smat.';
    T12G12=T_st12*BAst-G_st12*Smat.';
    T22G22=Tst_st22*BAst-G_st22*Smat.';
    X=[Tst_st11,T12G12;T_st21,T22G22]\[u01T12P0;u02T22P0];
    u1=X(1:length(face1)*3);
    P=X(length(face1)*3+1:end);
    %TGTG=T12G12/T22G22;
    %%
    %u1=(-TGTG*T_st21+Tst_st11)^-1*(u01T12P0-TGTG*u02T22P0);
    %u1=(-TGTG*T_st21+Tst_st11)\(u01T12P0-TGTG*u02T22P0);
    %P=T12G12\(u01T12P0-Tst_st11*u1);
end