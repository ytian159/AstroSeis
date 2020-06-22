function topos = smooth_topo_rand ( llon, llat, hh, lmax)
% theta = llat(:)*pi/180;
phi = llon(1,:)*pi/180;
theta = llat(:,1)*pi/180;
% from Yuan
% % load('topo_d120.mat');
% % maxtopol=10;
% % nl=maxtopol;
% % k1=(nl+1)^2;
% % topoc=coef_topo(1:k1);
% % nfold=1;
% % h1=topo_gens(theta,phi,topoc,maxtopol,nfold);
% % topos = reshape( h1, size(llon));
% % return

% dtheta= (llat(2,1)-llat(1,1))*pi/180;
% dphi = (llon(1,2)-llon(1,1))*pi/180;
dtheta = theta(2)-theta(1);
dphi = phi(2)-phi(1);
[pp,tt]= meshgrid(phi,theta);

k=0;
for ell = 0: lmax
    %     for m =-ell: ell
    for m =0: ell
        k = k + 1;
        basis(k).l = ell;
        basis(k).m = m ;
        ylm1 = ylms2(ell, m , theta, phi);
        ylm1 = ylm1(:);
        basis(k).coef = sum(hh(:).*conj(ylm1).*sin(tt(:)))*dtheta*dphi;
    end
end
% synthesis
nbasis = k;
topos = 0;
rng('default');
rng(3);
maxh=20;
coefrand=rand(nbasis,1);

%coefrand=coefrand.*exp(-(0:nbasis-1)/40).';
for j =1: nbasis
    basis(j).coef =(coefrand(j)-0.5)*maxh;
    ylm1 = ylms2(basis(j).l, basis(j).m , theta, phi);
    %ylm1=ylm1(:);
    fac =1; 
    if basis(j).m > 0, fac = 2; end 
    topos = topos + basis(j).coef*ylm1(:)*fac ;
end
topos = reshape(real(topos), size(llat));
return
