function [face,numface,ds,xs0,ys0,zs0,thetas,phis,height,V,Tri]=gen_mesh_ph_topo(Rm,n,nfold)

[V,Tri]=ParticleSampleSphere('N',n);
fv=struct('faces',Tri,'vertices',V); 
fv_new=SubdivideSphericalMesh(fv,1);
V=fv_new.vertices;
Tri=fv_new.faces;
V=V*Rm;
% figure; hold on;
% h= patch('faces',Tri,'vertices',V);
% alpha(h,.9);
% set(h,'EdgeColor','b','FaceColor',[1 1 1 ]*.5)
% axis equal vis3d
% view(3)
 %break;

%% get topo on patches 
for j =1: length(V)
    rtmp = norm(V(j,1:3));
    rtmp1 = norm(V(j,1:2));
    phi (j) = atan2( V(j,2), V(j,1));
    %     if phi(j) < 0, phi(j) = phi(j) + 360; end
    theta(j) = atan2(rtmp1,V(j,3));
end
rmars=Rm;

%% interp to get topography at the triangle vertex

% set topo
hh = get_phobos_topo (theta,phi );
hh=real(hh*nfold);
r0 = 10.993e3;
for j =1: length(V)
    r = (r0 + 1*hh(j));
    zj = r*cos(theta(j));
    xj = r*sin(theta(j))*cos(phi(j));
    yj = r*sin(theta(j))*sin(phi(j));
    V(j, 1:3) = [xj yj zj]; % vertex with topography
end
disp('... mesh on topo generated ');
%%
figure; hold on;
h= patch('faces',Tri,'vertices',V, 'FaceVertexCData', hh(:), 'FaceColor','interp');
colormap (jet);
alpha(h,.6);
set(h,'EdgeColor',[152 57 153]/255,'linewidth',.01);
% set(h,'EdgeColor','b','FaceColor',[1 1 1 ]*.5)
axis equal vis3d
view(3)
colorbar;
pause(.1);
%% set property of all triangles
numface = length(Tri);
np=numface;
sprintf('number of triangles = %d\n',numface)
pause(0.1);
xs0=zeros(np,1);
ys0=zeros(np,1);
zs0=zeros(np,1);
ds=zeros(np,1);
for j = 1:  numface
    ind = Tri(j,:);
    % 3 vertices of the triangle
    face(j).A = V(ind(1),:);
    face(j).B= V(ind(2),:);
    face(j).C = V(ind(3),:);
    % 3 sides
    a = norm( face(j).B- face(j).C);
    b = norm( face(j).C- face(j).A);
    c = norm( face(j).B- face(j).A);
    face(j).a = a;
    face(j).b = b;
    face(j).c = c;
    s = (a + b + c )/2;
    face(j).area = sqrt( s* (s- a)*(s- b)*(s- c));
    ds(j)=sqrt( s* (s- a)*(s- b)*(s- c));
    r = face(j).area/s;
    %     r = sqrt( (s-a)*(s- b)*(s-c)/s);
    face(j).r = r; % inradius
    % face normal
    vec_A_to_B = face(j).B- face(j).A;
    vec_A_to_B = vec_A_to_B/norm(vec_A_to_B);
    vec_A_to_C = face(j).C- face(j).A;
    vec_A_to_C = vec_A_to_C/norm(vec_A_to_C);
    vec_N = cross(vec_A_to_B, vec_A_to_C);
    vec_n = vec_N / norm(vec_N);
    if dot(vec_n, face(j).A) < 0 % set outward normal
        vec_n = - vec_n;
    end
    face(j).nvec = vec_n;
    % incenter location
    cost = dot ( vec_A_to_B, vec_A_to_C);
    sinhalft=sqrt( (1-cost)/2);
    dis = r /sinhalft;
    bivec = (vec_A_to_B + vec_A_to_C);
    bivec = bivec/norm(bivec);
    pos_ic = face(j).A + dis*bivec;
    face(j).ic = pos_ic;
    xs0(j)=pos_ic(1);
    ys0(j)=pos_ic(2);
    zs0(j)=pos_ic(3);
    %     plot3 ( [face(j).A(1) face(j).B(1) face(j).C(1) face(j).A(1)],...
    %         [face(j).A(2) face(j).B(2) face(j).C(2) face(j).A(2)],...
    %         [face(j).A(3) face(j).B(3) face(j).C(3) face(j).A(3)],'r-');
    %pend = pos_ic + scl/10* face(j).nvec  ;
    height(j)=norm(pos_ic)-Rm;
end
[phis,lat,rid]=cart2sph(xs0,ys0,zs0);
thetas=pi/2-lat;
%height(1:np)=0;
end