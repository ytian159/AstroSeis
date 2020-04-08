function ST = int_self_trac_tri_vec (vp,vs,rho,w,...
    xs,ys,zs,face,wi,ref,ridfac)
n1 = face.nvec(1);
n2 = face.nvec(2);
n3 = face.nvec(3);
%nxi = 300; 
% map xy(1:2, :) to real locations
for j =1: 3
    X(j,:) = face.A(j)*(1 - ref(1,:)-ref(2,:)) + ...
        face.B(j)* ref(1,:) + ...
        face.C(j)*ref(2,:);
end
xi = X(1,:);
yi = X(2,:);
zi = X(3,:);
rr = sqrt( (xi-face.ic(1)).^2 + (yi-face.ic(2)).^2+(zi-face.ic(3)).^2); 
ind2 = find(rr >= face.r*ridfac); 
xi = xi(ind2);
yi = yi(ind2); 
zi = zi(ind2); 

% figure; 
% plot3(xi,yi,zi,'.'); 

[T11, T12, T13,...
    T21,T22,T23,...
    T31, T32, T33] = green_traction_tensor (vp,vs,rho,w,...
    xi,yi,zi,xs,ys,zs,n1,n2,n3);
fac = 2; 
st11 = sum(T11.*wi) *face.area*fac;
st12 = sum(T12.*wi) *face.area*fac;
st13 = sum(T13.*wi) *face.area*fac;
st21 = sum(T21.*wi) *face.area*fac;
st22 = sum(T22.*wi) *face.area*fac;
st23 = sum(T23.*wi) *face.area*fac;
st31 = sum(T31.*wi) *face.area*fac;
st32 = sum(T32.*wi) *face.area*fac;
st33 = sum(T33.*wi) *face.area*fac;
ST =[st11 st12 st13
    st21 st22 st23
    st31 st32 st33];
 ST = ST -1/2 * eye(3); 
return
end
