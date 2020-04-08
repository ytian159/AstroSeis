function [P]=asso_P_theta0(l,m,theta)
P=1/2^m*factorial(l+m)/factorial(l-m)*theta.^m;

end