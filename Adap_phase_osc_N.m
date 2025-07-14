function [dvar] = Adap_phase_osc_N(var,par,N)
%Adap_phase_osc_N: the ODe function for adaptive phase oscillator system 
%   inputs:
%   var: system variables vector of length N+1
%   par: system parameters
%   N: number of oscillators in the network 
%   output: 
%   dvar: the time dervative for the variables


% variables
phi = var(1:N);
mu = var(N+1);

% parameters
ome = par(1:N);
kappa = par(N+1)./N;
eta = par(N+2);
alpha = par(N+3);
EPSS = par(N+4);

X = (1./N).*sum( sin(phi+alpha) );

% producing the coupling terms by repeating the phi vector into a metrix
% substract the transpose and then taking sin and summing up the colomns

coupling_vector = sum(sin(repmat(phi,1,N) - repmat(phi,1,N)'));

dvar(1:N) = (ome + mu) - sin(phi) + kappa.*(coupling_vector');  
dvar(N+1) = EPSS*(-mu + eta*(1-X));
dvar = dvar';

end