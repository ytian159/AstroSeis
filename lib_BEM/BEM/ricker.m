function [w,t] = ricker(f0,ts,dt,nt)
% [w,t] = ricker(f0,ts,dt,nt)
% W(t) = ( 1-2at^2) exp ( -at^2), a=(pi*t)^2

% f0: central frequency in Hz
% ts: shift in sec 
% dt: time sampling inteval in sec.
% nt: number of time samples 

a = (pi*pi)*(f0*f0);
t = [0:nt-1]*dt;
tp = t - ts;
arg = a*tp.^2;
w = ( 1. - 2*a*tp.*tp).*exp(-arg);
