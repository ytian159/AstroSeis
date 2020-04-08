function [Smat]=Smat_func(face)
ngd=length(face);
Smat=zeros(ngd,3*ngd);
for i=1:ngd
    Smat(i,i)=face(i).nvec(1);
    Smat(i,i+ngd)=face(i).nvec(2);
    Smat(i,i+2*ngd)=face(i).nvec(3);
    
end

end