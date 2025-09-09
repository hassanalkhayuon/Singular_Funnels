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

%% settings 
Eps_res = 1000;
Del_res = 1000;

Eps_arr = linspace(0.02,0.08,Eps_res);
Del_arr = linspace(0.6,1.6,Del_res);

%% Parameter scan 
Om1byOm2 = NaN(Del_res,Eps_res);

% Loop over epsilon values
for ind_Eps = 1:Eps_res
    Eps = Eps_arr(ind_Eps);
    parfor ind_Del = 1:Del_res
        Del = Del_arr(ind_Del);
        Om1byOm2(ind_Del,ind_Eps) = Ome1byOme2_fun(Del,Eps);
    end
    disp(Eps_res - ind_Eps)
end

%% 1:1 locking curve
res_1to1 = 100;
Eps_1to1 = linspace(0.02,0.08,res_1to1);
Del_1to1 = NaN(size(Eps_1to1));
Del_init = 1;
for ind_1to1 = 1:res_1to1
    Eps = Eps_1to1(ind_1to1);
    gg = @(DEL)(Ome1byOme2_fun(DEL,Eps)-1.01);
    Del_init = fzero(gg,Del_init);
    Del_1to1(ind_1to1) = Del_init;
    disp([res_1to1 - ind_1to1,Ome1byOme2_fun(Del_init,Eps)]);
end
figure(1)
hold on 
plot(Eps_1to1,Del_1to1,'--w')

%%
figure(10)
clf
Scan_plot = pcolor(Eps_arr,Del_arr,Om1byOm2);
Scan_plot.LineStyle = "none";
clim([min(min(Om1byOm2)), max(max(Om1byOm2))])
set(gca,'FontSize',15)
xlabel('$\varepsilon$')
ylabel('$\Delta_{\omega}$','Rotation',0)


%%
function output = Ome1byOme2_fun(Del,Eps)

N = 2;
eta = 10; % adaptive parameters
alpha = pi/2; % phase shift
kappa = 1;
ome2 = -4;

opts = odeset('RelTol',1e-5,'AbsTol',1e-5);

tend = 100 / Eps;

ome = [ome2 + Del; ome2];
par = [ome; kappa; eta; alpha; Eps];

odefun = @(t, var) Adap_phase_osc_N(var, par, N);

[~, var_temp] = ode15s(odefun, [0 tend], [pi; pi/2; 15], opts);
[~, var] = ode15s(odefun, [0 tend], var_temp(end,:), opts);

phi1 = var(:,1);
phi2 = var(:,2);

Om1= (phi1(1) - phi1(end)) / tend;
Om2= (phi2(1) - phi2(end)) / tend;

output = Om1/Om2;
end
