function [SBJ]=SphericalBesselJ(l,kr)

SBJ=besselj(l+1/2,kr).*sqrt(pi./(2.*kr));
end