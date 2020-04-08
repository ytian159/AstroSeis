function [GRT]=G_singular_value(vp,vs,rho,w,...
    R,n1,n2,n3)
kp = w/vp; 
ks = w/vs;

value1=(1/4).*R.^(-1).*((sqrt(-1)*2).*ks.*R+exp(1).^(sqrt(-1).*kp.*R).*( ...
  1+(sqrt(-1)*(-1)).*kp.*R)+exp(1).^(sqrt(-1).*ks.*R).*((-1)+(sqrt( ...
  -1)*(-1)).*ks.*R))./(rho.*w^2);

value2=(1/2).*R.^(-1).*((-1).*exp(1).^(sqrt(-1).*kp.*R)+sqrt(-1).*kp.*R+ ...
  exp(1).^(sqrt(-1).*ks.*R).*(1+(sqrt(-1)*(-1)).*ks.*R))./(rho.*w^2);
Galy=[value1,0,0;0,value1,0;0,0,value2];
np=[n1,n2,n3];
phi_rot=atan2(np(2),np(1));
theta_rot=atan2(sqrt(np(1)^2+np(2)^2),np(3));
R_mat=rotz(phi_rot*180/pi)*roty(theta_rot*180/pi); % This matrix is to Galy rotate the matrix from local coordinate to cartisan coor
%%

%n1*nx+n2*ny+n3*nz=0
% R1=[n1/sqrt(n1^2+n2^2),-n2/sqrt(n1^2+n2^2),0;n2/sqrt(n1^2+n2^2),n1/sqrt(n1^2+n2^2),0;0,0,1];
% R2=[n3/sqrt(n1^2+n3^2),0,n1/sqrt(n1^2+n3^2);0,1,0;-n1/sqrt(n1^2+n3^2),0,n3/sqrt(n1^2+n3^2)];
% HRT=H0*R1*R2;
GRT=R_mat*Galy*R_mat.';

end