function ST = int_self_A_tri_vec (vp,w,...
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

Gpn = greens_func_deri_p (vp,w,...
    xi,yi,zi,xs,ys,zs,n1,n2,n3);
fac = 2; 
ST = sum(Gpn.*wi) .*face.area.*fac;
ST = ST -1/2; 
%ST=-1/2;

return
end
