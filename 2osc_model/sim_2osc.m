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
initcond = [1.2976; 5.06437; -0.137124];
% initcond = [1.2461; 3.2; 3.8];


% Loop over epsilon values
Epss = 0.08;
tend = 100 / Epss;

Del = 1.5;
ome = [ome2 + Del; ome2];

par = [ome; kappa; eta; alpha; Epss];

odefun = @(t, var) Adap_phase_osc_N(var, par, N);
[~, var] = ode15s(odefun, [0 tend], initcond, opts);
% [~, var] = ode15s(odefun, [0 tend], var_temp(end,:), opts);

mu = var(:,3);
phi1 = mod(var(:,1),2*pi);
phi2 = mod(var(:,2),2*pi);

%%
figure
plot3(mu,phi1,phi2,'.k')
xlim([0,2*pi])
ylim([0,2*pi])
%%
% Scan_plot = pcolor(Eps_arr,Del_arr,XXX);
% Scan_plot.LineStyle = "none";
% set(gca,'FontSize',15)
% xlabel('$\varepsilon$')
% ylabel('$\Delta_{\omega}$','Rotation',0)