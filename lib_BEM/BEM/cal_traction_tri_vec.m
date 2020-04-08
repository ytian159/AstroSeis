function [TRAC2]=cal_traction_tri_vec(face,numface,w,lamda,mu,rho,Q)
%% Vectorize the traction matrix calculation
% Yuan Tian 04/01/2017 at UH
vp0=sqrt((lamda+2*mu)/rho);
vs0=sqrt(mu/rho);
Qp=2.5*Q;
vp=vp0/(1+complex(0,1)*0.5/Qp);
vs=vs0/(1+complex(0,1)*0.5/Q);
%% Creating sub mesh
 ngd=length(face);
 nint = 10; % This number can only be 10, 20, 40, 50
 %nsub=nint^2/4;
 [xi,yi,zi,n1,n2,n3,wi]=subgrid_gen(face(1),nint);
 nsub=length(xi);
 xia=zeros(numface,nsub);
  yia=zeros(numface,nsub);
zia=zeros(numface,nsub); 
wia=zeros(numface,nsub); 
 n1a=zeros(numface,1);
 n2a=zeros(numface,1);
 n3a=zeros(numface,1);
 area=zeros(numface,1);
 for i=1:ngd
     [xi,yi,zi,n1,n2,n3,wi]=subgrid_gen(face(i),nint);
     xia(i,:)=xi;
     yia(i,:)=yi;
     zia(i,:)=zi;
     n1a(i)=n1;
     n2a(i)=n2;
     n3a(i)=n3;
     wia(i,:)=wi;
     area(i)=face(i).area;
 end
%% Self element integral
STs=zeros(3,3,numface);
nxi=300;
ridfac=1;
xi = linspace(0,1, nxi); 
dx = xi(2)-xi(1) ;
etai = linspace(0,1,nxi);
[xx,nn]=meshgrid(xi,etai);
xi1 = xx(:);
eta1=nn(:);
ind = find( xi1 + eta1 <= 1); 
xi1 = xi1(ind);
eta1 = eta1(ind); 
refs(1,:) = xi1; 
refs(2,:) = eta1; 
wis = dx * dx; 
parfor i=1:numface
    ps=face(i).ic;
    xs=ps(1);
    ys=ps(2);
    zs=ps(3);
    STs(:,:,i)= int_self_trac_tri_vec (vp,vs,rho,w,...
        xs,ys,zs,face(i),wis,refs,ridfac)+eye(3);
end
%% Overall integral and constructing matrix
%starray=zeros(numface,numface,3,3);
TRAC2=zeros(3*numface,3*numface);
for j=1:numface
     ps=face(j).ic;
 xs=ps(1);
 ys=ps(2);
 zs=ps(3);
    ST0=trac_tri_int_vec(vp,vs,rho,w,xia,yia,zia,wia,xs,ys,zs,area,n1a,n2a,n3a,nsub,STs(:,:,j),j);
    TRAC2(j,:)=ST0(:,1);
    TRAC2(j+numface,:)=ST0(:,2);
    TRAC2(j+2*numface,:)=ST0(:,3);
end

end