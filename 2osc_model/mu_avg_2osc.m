% This code is created by Hassan Alkhayuon to
% calculate the avreaging system for mu in the 2 pahse oscilator system 
% A Research project with Serhiy Yanchuk and Hildeberto Jardón-Kojakhmetov
% and Sebastian Wieczorek 


function [g_bar] = mu_avg_2osc(mu,par)
% This is a avrage dynamics function for the avrage system for 2 coupled
% adaptive phase oscillators

% System parameters

ome = par(1:2);
kappa = par(3);
eta = par(4);
alpha = par(5);

par = [ome;  kappa; eta; alpha; mu];

% numarical parameters 
Ttrans = 100;
Tend = 500;
opts = odeset('RelTol',1e-5,'AbsTol',1e-5);

% initial condion for fast subsystem 
phi1_init = rand*2*pi;
phi2_init = rand*2*pi;
initcond = [phi1_init; phi2_init];

% ode function
odefun = @(t,var)Adap_phase_osc_2fast(var,par);

% first integration to remove the transiant 
[~,phi_temp] = ode45(odefun, [0, Ttrans], initcond, opts);

% update the initial condition 
initcond = mod(phi_temp(end,:),2*pi);

% second integration 
[t,phi] = ode45(odefun, 0:0.1:Tend, initcond, opts);

X = (1/2) * sum( sin(phi + alpha), 2 ); % Sum across columns (dimension 2)
gt = -mu + eta * (1 - X);

% Avraging 
g_bar = trapz(t,gt)./Tend; 

end


