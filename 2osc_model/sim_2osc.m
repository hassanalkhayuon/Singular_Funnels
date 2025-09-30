% This code is created by Hassan Alkhayuon to
% a parameter plane for the 2 coupled adaptive phase
% oscillators to show syncrnisation vs nonsyncrnisation regions
% A Research project with Serhiy Yanchuk, Hildeberto Jardón-Kojakhmetov
% and Sebastian Wieczorek

clear
warning off
addpath("../")

%% Parameters
N = 2;
eta = 10; % adaptive parameters
alpha = pi/2; % phase shift
kappa = 1;
ome2 = -4;

opts = odeset('RelTol',1e-5,'AbsTol',1e-5);

%% Simulations


% initcond = [pi; pi/2; 15];
initcond = [0; 0; 10.0302];
% initcond = [1.2461; 3.2; 3.8];


% Loop over epsilon values
Epss = 0.1;
tend = 1.1;

Del = 1.5;
ome = [ome2 + Del; ome2];

par = [ome; kappa; eta; alpha; Epss];

odefun = @(t, var) Adap_phase_osc_N(var, par, N);
[t, var] = ode15s(odefun, [0 tend], initcond, opts);
% [~, var] = ode15s(odefun, [0 tend], var_temp(end,:), opts);

mu = var(:,3);
phi1 = mod(var(:,1),2*pi);
phi2 = mod(var(:,2),2*pi);

%%
figure
plot(mu,var(:,2), '-k')
% xlim([0,2*pi])
% ylim([0,2*pi])
%%
% Scan_plot = pcolor(Eps_arr,Del_arr,XXX);
% Scan_plot.LineStyle = "none";
% set(gca,'FontSize',15)
% xlabel('$\varepsilon$')
% ylabel('$\Delta_{\omega}$','Rotation',0)