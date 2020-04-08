function [TRAC2]=cal_A_st(face1,face2,w,lamda,mu,rho,Q,nint,nxi)
%% Vectorize the traction matrix calculation
% Yuan Tian 04/01/2017 at UH
vp0=sqrt((lamda+2*mu)/rho);
vs0=sqrt(mu/rho);
Qp=2.5*Q;
vp=vp0/(1+complex(0,1)*0.5/Qp);
vs=vs0/(1+complex(0,1)*0.5/Q);
%% Creating sub mesh
ngd1=length(face1);
ngd2=length(face2);
%nint = 10; % This number can only be 10, 20, 40, 50
%nsub=nint^2/4;
[xi,yi,zi,n1,n2,n3,wi]=subgrid_gen(face2(1),nint);
nsub=length(xi);
xia=zeros(ngd2,nsub);
yia=zeros(ngd2,nsub);
zia=zeros(ngd2,nsub);
wia=zeros(ngd2,nsub);
n1a=zeros(ngd2,1);
n2a=zeros(ngd2,1);
n3a=zeros(ngd2,1);
area=zeros(ngd2,1);
for i=1:ngd2
    [xi,yi,zi,n1,n2,n3,wi]=subgrid_gen(face2(i),nint);
    xia(i,:)=xi;
    yia(i,:)=yi;
    zia(i,:)=zi;
    n1a(i)=n1;
    n2a(i)=n2;
    n3a(i)=n3;
    wia(i,:)=wi;
    area(i)=face2(i).area;
end
%% Self element integral
if ngd1==ngd2
    STs=zeros(ngd2,1);
    %nxi=200;
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
    parfor i=1:ngd2
        ps=face2(i).ic;
        xs=ps(1);
        ys=ps(2);
        zs=ps(3);
        STs(i)= int_self_A_tri_vec (vp,w,...
            xs,ys,zs,face2(i),wis,refs,ridfac)+1;
    end
    
%% Overall integral and constructing matrix

    TRAC2=zeros(ngd1,ngd2);
    for j=1:ngd1
        ps=face2(j).ic;
        xs=ps(1);
        ys=ps(2);
        zs=ps(3);
        ST0=A_tri_int_vec(vp,w,xia,yia,zia,wia,xs,ys,zs,area,n1a,n2a,n3a,nsub,STs(j),j);
        TRAC2(j,:)=ST0;
    end
else
   TRAC2=zeros(ngd1,ngd2);
    for j=1:ngd1
        ps=face1(j).ic;
        xs=ps(1);
        ys=ps(2);
        zs=ps(3);
        ST0=A_tri_int_vec_diflay(vp,w,xia,yia,zia,wia,xs,ys,zs,area,n1a,n2a,n3a,nsub);
        TRAC2(j,:)=ST0;
    end 
    
end
    
    
    
end