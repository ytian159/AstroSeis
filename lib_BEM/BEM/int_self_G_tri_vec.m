function ST = int_self_G_tri_vec (vp,vs,rho,w,...
    xs,ys,zs,face,wi,ref,ridfac)
n1 = face.nvec(1);
n2 = face.nvec(2);
n3 = face.nvec(3);
%nxi = 300; 
% map xy(1:2, :) to real locations
np=length(ref);
X=zeros(3,np);
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

[G11, G12, G13,...
    G21,G22,G23,...
    G31, G32, G33] = Bmat_func (vp,vs,rho,w,...
    xi,yi,zi,xs,ys,zs);
fac = 2; 
st11 = sum(G11.*wi) *face.area*fac;
st12 = sum(G12.*wi) *face.area*fac;
st13 = sum(G13.*wi) *face.area*fac;
st21 = sum(G21.*wi) *face.area*fac;
st22 = sum(G22.*wi) *face.area*fac;
st23 = sum(G23.*wi) *face.area*fac;
st31 = sum(G31.*wi) *face.area*fac;
st32 = sum(G32.*wi) *face.area*fac;
st33 = sum(G33.*wi) *face.area*fac;
ST =[st11 st12 st13
    st21 st22 st23
    st31 st32 st33];
 ST=ST+G_singular_value(vp,vs,rho,w,...
   face.r,n1,n2,n3);
%  ST=G_singular_value(vp,vs,rho,w,...
%      face.r,n1,n2,n3);
 return
end
