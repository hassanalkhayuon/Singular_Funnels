% This code is created by Hassan Alkhayuon to
% the eps-eta curve for homoclinic cycle whihc forms the boundary of 
% singuler basin region.
% This is a Research project with Serhiy Yanchuk 
% and Hildeberto Jardón-Kojakhmetov

clear
warning off
addpath('..\')

% Colour-blind friendly colour
[red,yellow,green,blue ] = ...
    Colour_blind_friendly_colours();

%% Parameters
options = optimset('TolX',1e-7);

eta = 10; % adaptive parameters
ome = -4;
alpha0 = 2.62281789592832; % phase shift.

EPS_arr = linspace(0.1,0.001,10);
for ind = 1:length(EPS_arr)
    EPS = EPS_arr(ind);
    gap_fun = @(alpha)homoclinic_gap_e2(alpha,EPS,ome,eta);
    alpha0 = fzero(gap_fun,alpha0, options);
    alpha_arr(ind) = alpha0;
    disp([alpha_arr(ind),ind])
end








%% 
function [output] = homoclinic_gap_e2(alpha,EPS,ome,eta)
par = [ome, eta, alpha, EPS];


% equilibria
phi_sym = sym('phi_sym');

eq_expr = -ome + sin(phi_sym) == eta*(1- sin(phi_sym+alpha));

phi_eq = solve(eq_expr, phi_sym,Real=true);

phi_e1 = mod(eval(phi_eq(2)),2*pi);
phi_e2 = mod(eval(phi_eq(1)),2*pi);

mu_e1 = sin(phi_e1) - ome;

mu_e2 = sin(phi_e2) - ome;


e1 = [phi_e1, mu_e1];

opts = odeset('RelTol',1e-10,'AbsTol',1e-10);

odefun = @(t,var)Adaptive_phase_ode(var,par);


% figure(4);
% cla
% hold on
% axis([1 5 -2*pi 2*pi])
% 
initcond_p = [phi_e2 mu_e2] + 0.001.*[0 1];

[~, var] = ode15s(odefun,[1000 0],initcond_p,opts);
    
% phi = var(:,1);
% mu = var(:,2);

% plot(mu, phi,'-','Color','b','LineWidth',2)

if abs(var(end,2) - mu_e1) < 0.01
    output = 1;
else
    output = -1;
end

end

%% ODE function
function [dvar] = Adaptive_phase_ode(var,par)

% model parameters

I       = par(1);
eta     = par(2);
alpha   = par(3);
EPS     = par(4);


% variables
phi = var(1);
mu = var(2);

% diffrential equations
dphi = I + mu - sin(phi);
dmu  = EPS*(-mu + eta*(1 - sin(phi + alpha)));

dvar = [dphi; dmu];
end