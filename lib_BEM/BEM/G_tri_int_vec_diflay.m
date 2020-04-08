function [ST,Gdum]=G_tri_int_vec_diflay(vp,vs,rho,w,...
    xi,yi,zi,wi,xs,ys,zs,area)

[T11, T12, T13,...
    T21,T22,T23,...
    T31, T32, T33] = Bmat_func (vp,vs,rho,w,...
    xi,yi,zi,xs,ys,zs);
fac = 2; 
% st11 = sum(sum(T11.*wi,2) .*area)*fac;
% st12 = sum(sum(T12.*wi,2) .*area)*fac;
% st13 = sum(sum(T13.*wi,2) .*area)*fac;
% st21 = sum(sum(T21.*wi,2) .*area)*fac;
% st22 = sum(sum(T22.*wi,2) .*area)*fac;
% st23 = sum(sum(T23.*wi,2) .*area)*fac;
% st31 = sum(sum(T31.*wi,2) .*area)*fac;
% st32 = sum(sum(T32.*wi,2) .*area)*fac;
% st33 = sum(sum(T33.*wi,2) .*area)*fac;
st11 = sum(T11.*wi,2) .*area*fac;
st12 = sum(T12.*wi,2) .*area*fac;
st13 = sum(T13.*wi,2) .*area*fac;
st21 = sum(T21.*wi,2) .*area*fac;
st22 = sum(T22.*wi,2) .*area*fac;
st23 = sum(T23.*wi,2) .*area*fac;
st31 = sum(T31.*wi,2) .*area*fac;
st32 = sum(T32.*wi,2) .*area*fac;
st33 = sum(T33.*wi,2) .*area*fac;


% STs = int_self_green_traction_tensor_triangle (vp,vs,rho,w,...
%                 xs,ys,zs,face);
% STs=STs+eye(3
ST =[st11, st12, st13;
    st21, st22, st23;
    st31, st32, st33];
Gdum=[st11(1),st12(1),st13(1);st21(1),st22(1),st23(1);st31(1),st32(1),st33(1)];
% ST =[st11 st21 st31
%     st12 st22 st32
%     st13 st23 st33];
% ST(:,1,1)=st11;
% ST(:,1,2)=st12;
% ST(:,1,3)=st13;
% 
% ST(:,2,1)=st21;
% ST(:,2,2)=st22;
% ST(:,2,3)=st23;
% 
% ST(:,3,1)=st31;
% ST(:,3,2)=st32;
% ST(:,3,3)=st33;
end