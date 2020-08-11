function [T11, T12, T13,T21,T22,T23,T31, T32, T33]= green_traction_tensor_simp (vp,vs,rho,w,...
    x1,x2,x3,xs1,xs2,xs3,sn1,sn2,sn3)
% Tij: traction vector i-th component; for the single-force 
% direction along 'j-th' direction; 
% vp, vs, rho: P-wave, s-wave velocity and density 
% w: angular frequency = 2*pi*f 
% (x1,x2,x3) field point 
% (sn1,sn2,sn3) surface element normal direction at the field point 
% (xs,ys,zs) source point 
im = complex(0,1); 
R = sqrt( (x1-xs1).^2 + (x2-xs2).^2 + (x3-xs3).^2); 
gm1 = (x1-xs1)./R; 
gm2 = (x2-xs2)./R; 
gm3 = (x3-xs3)./R; 
kp = w/vp; 
ks = w/vs;
gamma=vs/vp;
xP=kp.*R;
xS=ks.*R;
B=(3./xS.^2-3.*im./xS-1).*exp(xS.*im)-gamma.^2.*(...
    3./xP.^2-3.*im./xP-1).*exp(xP.*im);
C=(-15./xS.^2+15.*im./xS+6-im*xS).*exp(xS.*im)-gamma.^2.*(...
    -15./xP.^2+15.*im./xP+6-im*xP).*exp(xP.*im);
D=(im.*xS-1).*exp(xS.*im)+2*B;
E=(1-2*gamma^2)*(im*xP-1).*exp(xP.*im)+2*B;
coef=1./(4*pi*R.^2);
gm11C=2*C.*gm1.^2;
gm22C=2*C.*gm2.^2;
gm33C=2*C.*gm3.^2;
gm123C=2*C.*gm1.*gm2.*gm3;
gm122C=gm22C.*gm1;
gm133C=gm33C.*gm1;
gm112C=gm11C.*gm2;
gm113C=gm11C.*gm3;
gm223C=gm22C.*gm3;
gm233C=gm33C.*gm2;
gm1D=(gm1).*D;
gm2D=(gm2).*D;
gm3D=(gm3).*D;
gm1E=(gm1).*E;
gm2E=(gm2).*E;
gm3E=(gm3).*E;
%% Tijk=coef*[2*C*gmi*gmk*gmj+(Dtik*gmj+Dtjk*gmi)D+Dtij*gmk*E]
T111=(gm11C.*gm1+2.*gm1D+gm1E);
T121=(gm112C+gm2D);
T131=(gm113C+gm3D);

T211=(gm112C+gm2D);
T221=(gm122C+gm1E);
T231=(gm123C);

T311=(gm113C+gm3D);
T321=(gm123C);
T331=(gm133C+gm1E);

T112=(gm112C+gm2E);
T122=(gm122C+gm1D);
T132=(gm123C);

T212=(gm122C+gm1D);
T222=(gm22C.*gm2+2.*gm2D+gm2E);
T232=(gm223C+gm3D);

T312=(gm123C);
T322=(gm223C+gm3D);
T332=(gm233C+gm2E);

T113=(gm113C+gm3E);
T123=(gm123C);
T133=(gm133C+gm1D);

T213=(gm123C);
T223=(gm223C+gm3E);
T233=(gm233C+gm2D);

T313=(gm133C+gm1D);
T323=(gm233C+gm2D);
T333=(gm33C.*gm3+2.*gm3D+gm3E);

%% Tik=Tijk*nj
T11=coef.*(T111.*sn1+T121.*sn2+T131.*sn3);
T12=coef.*(T112.*sn1+T122.*sn2+T132.*sn3);
T13=coef.*(T113.*sn1+T123.*sn2+T133.*sn3);

T21=coef.*(T211.*sn1+T221.*sn2+T231.*sn3);
T22=coef.*(T212.*sn1+T222.*sn2+T232.*sn3);
T23=coef.*(T213.*sn1+T223.*sn2+T233.*sn3);

T31=coef.*(T311.*sn1+T321.*sn2+T331.*sn3);
T32=coef.*(T312.*sn1+T322.*sn2+T332.*sn3);
T33=coef.*(T313.*sn1+T323.*sn2+T333.*sn3);


end