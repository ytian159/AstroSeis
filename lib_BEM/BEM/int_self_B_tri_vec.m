function ST = int_self_B_tri_vec (vp,w,...
    xs,ys,zs,face,wi,ref,ridfac)

%nxi = 300; 
% map xy(1:2, :) to real locations
kp = w/vp; 

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
%%
Gpn = greens_func_p (vp,w,...
    xi,yi,zi,xs,ys,zs);
fac = 2; 
ST = sum(Gpn.*wi) .*face.area.*fac;
R=face.r;
B_singular=(sqrt(-1)*(-1/2)).*((-1)+exp(1).^(sqrt(-1).*kp.*(R.^2).^(1/2))).* ...
  kp.^(-1);
ST=ST+B_singular;
% rest_area=sum(wi*length(xi)) .*face.area.*fac
% circle_area=pi*face.r^2
% total_area=face.area
%numerical_area=
%ST=B_singular;

return
end