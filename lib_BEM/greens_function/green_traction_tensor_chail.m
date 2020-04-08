function [T]= green_traction_tensor_chail (vp,vs,rho,w,...
    x1,x2,x3,xs1,xs2,xs3,sn1,sn2,sn3)
im = complex(0,1); 
ngd=length(x1);
R = sqrt( (x1-xs1).^2 + (x2-xs2).^2 + (x3-xs3).^2); 
gm(:,1) = (x1-xs1)./R; 
gm(:,2) = (x2-xs2)./R; 
gm(:,3) = (x3-xs3)./R; 
kp = w/vp; 
ks = w/vs;
gamma=vs/vp;
xP=kp.*R;
xS=ks.*R;
delta=eye(3);
B=(3./xS.^2-3.*im./xS-1).*exp(xS.*im)-gamma.^2.*(...
    3./xP.^2-3.*im./xP-1).*exp(xP.*im);
C=(-15./xS.^2+15.*im./xS+6-im*xS).*exp(xS.*im)-gamma.^2.*(...
    -15./xP.^2+15.*im./xP+6-im*xP).*exp(xP.*im);
D=(im.*xS-1).*exp(xS.*im)+2*B;
E=(1-2*gamma^2)*(im*xP-1).*exp(xP.*im)+2*B;
coef=1./(4*pi*R.^2*rho*w^2);
Tijk=zeros(ngd,3);
T=zeros(ngd,3,3);
for k=1:3
    for i=1:3
        for j=1:3
            Tijk(:,j)=coef*(2*C.*gm(:,i).*gm(:,k).*gm(:,j)+...
                (delta(i,k).*gm(:,j)+delta(j,k).*gm(:,i)).*D...
                +delta(i,j).*gm(:,k).*E);
        end
        T(:,k,i)=Tijk(:,1).*sn1+Tijk(:,2).*sn2+Tijk(:,3).*sn3;
    end
end

end