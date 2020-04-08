function [xi,yi,zi,n1,n2,n3,wi]=subgrid_gen(face,n)

n1 = face.nvec(1);
n2 = face.nvec(2);
n3 = face.nvec(3);
%n = 10; 
numnodes = rule_full_size ( n )  ; 
vert1 =[0 0 ]; 
vert2=[1 0]; 
vert3=[ 0 1]; 
[ ref, wi ] = triasymq ( n, vert1, vert2, vert3, numnodes ); 
np=length(ref);
% map xy(1:2, :) to real locations
X=zeros(3,np);
for j =1: 3
    X(j,:) = face.A(j)*(1 - ref(1,:)-ref(2,:)) + ...
        face.B(j)* ref(1,:) + ...
        face.C(j)*ref(2,:);
end
xi = X(1,:);
yi = X(2,:);
zi = X(3,:);

end