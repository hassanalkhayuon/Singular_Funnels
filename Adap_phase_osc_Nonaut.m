function [dvar] = Adap_phase_osc_Nonaut(t,var,par)
%Adap_phase_osc_N: the ODE function for nonautonomous adaptive phase 
%  oscillator system 
%   inputs:
%   var: system variables vector of length N+1
%   par: system parameters
%   N: number of oscillators in the network 1 in this case 
%   output: 
%   dvar: the time dervative for the variables

N = 1;

% variables
phi = var(1);
mu = var(2);

% parameters
ome_start = par(1);
ome_del = par(2);
rate = par(3);
eta = par(4);
alpha = par(5);
EPSS = par(6);

ome_t = ome_start + (ome_del/2)*(tanh(rate*t) + 1);

dvar(1:N) = (ome_t + mu) - sin(phi);  
dvar(N+1) = EPSS*(-mu + eta*(1-sin(phi+alpha)));
dvar = dvar';

end