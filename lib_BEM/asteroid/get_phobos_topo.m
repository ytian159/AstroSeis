function topo = get_phobos_topo (theta,phi )
%% phobos topography
% based on the paper
% "Phobos' shape and topography models" by Willner et al., PSS paper 2014
% yingcai zheng, March 28, 2017
% yc.zheng@gmail.com
% %
% clear;close all; clc;
% set(0,'DefaultAxesFontSize',16);
% set(0,'defaultlinelinewidth',1)

r0 = 10.993e3; % mean phobos radius

aa= load('TABLEA1.DAT');

ll = aa(:,1);
mm = aa(:,2);
amn = aa(:,3);
bmn = aa(:,5);

% th = linspace(0, pi, 81);
% ph = linspace(-pi, pi, 111);
% [pp, tt] = meshgrid(ph, th);
pp = phi ;
tt = theta;

vv = zeros(size(pp));

for j = 1: numel(ll)
%     disp(j)
    fac =1;
    if mod(mm(j),2)==1
        fac=-1;
    end
    %     ctopo = ylm(ll(j), mm(j),th, ph)*fac;
    ctopo = fac * ylm1(ll(j),mm(j),tt(:),pp(:));
%     ctopo = reshape(ctopo, size(vv));
    vv = vv + amn(j)*real(ctopo) + bmn(j)*imag(ctopo);
end

topo = vv; 
return;

zz =  cos(tt);
xx =  sin(tt) .* cos(pp);
yy =  sin(tt).* sin(pp);

rr = r0 + vv*2.0;
zz1 =  rr.*cos(tt);
xx1 =  rr.*sin(tt) .* cos(pp);
yy1 =  rr.* sin(tt).* sin(pp);

surf(xx1,yy1,zz1,vv);
colormap jet;
shading interp;
% imagesc(ph, th, yy);
axis equal; hold on;
[mrow, ncol]= size(xx1);
for j =1: 5: ncol
    plot3(xx1(:,j), yy1(:,j), zz1(:,j),'color',[ 1 1 1]*.3);
end
for j =1: 5: mrow
    plot3(xx1(j,:), yy1(j,:), zz1(j,:),'color',[ 1 1 1]*.3);
end

figure;
vv = vv -r0 ;
pcolor(pp*180/pi, tt*180/pi, vv);
colormap jet;
hc=colorbar;
shading interp;
set(gca,'ydir','rev');
caxis([-3000 3000]);
ylabel(hc, 'Surface Topography (m)');

return
