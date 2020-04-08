function [ Green ] = greens_functionQ( w,rho,mu,lamda,x0,x,y,z,Q )
%UNTITLED7 ????????????
%   ????????
%x0 need to be dim 3 vector;
im_u=complex(0,1);
    Qs=Q;
    Qp=Qs*9/4;
    kp=w/sqrt((2*mu+lamda)/rho)*(1+im_u/(2*Qp));
    ks=w/sqrt((mu)/rho)*(1+im_u/(2*Qs));
dim=3;
%z=fsurf(x,y); %surfave function

r=sqrt((x0(1)-x)^2+(x0(2)-y)^2+(x0(3)-z)^2);
xr(1)=x-x0(1);
xr(2)=y-x0(2);
xr(3)=z-x0(3);
ts_prar=exp(sqrt(-1)*ks*r)/(4*pi*r^3);
tp_prar=exp(sqrt(-1)*kp*r)/(4*pi*r^3);
gamma=xr/r;
green_coef=1/(rho*w^2);
sigma=eye(dim);
gs2=zeros(3,3);
gp2=zeros(3,3);
Green=zeros(3,3);
for i=1:dim    
    for j=1:dim
        gs=ts_prar*r^2;
        gs2(i,j)=ts_prar*(((3-ks^2*r^2)*gamma(i)*gamma(j)-sigma(i,j))...
            +im_u*ks*r*(sigma(i,j)-3*gamma(i)*gamma(j)));
        gp2(i,j)=tp_prar*(((3-kp^2*r^2)*gamma(i)*gamma(j)-sigma(i,j))...
            +im_u*kp*r*(sigma(i,j)-3*gamma(i)*gamma(j)));
        Green(i,j)=green_coef*(gs*ks^2*sigma(i,j)+gs2(i,j)-gp2(i,j));
    end
end



end

