function [P]=asso_P_theta0_coef(l,m)
mabs=abs(m);
P=(-1/2)^mabs*((2*l+1)/(4*pi)*factorial(l+mabs)/factorial(l-mabs))^0.5/factorial(mabs);
if m<0
    P=conj(P)*(-1)^mabs;
end
end