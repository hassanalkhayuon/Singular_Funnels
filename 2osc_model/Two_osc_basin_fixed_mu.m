% This code is created by Hassan Alkhayuon to
% compute basin of attraction for 2 coupled adaptive phase oscillators,
% with fixed mu.
% A Research project with Serhiy Yanchuk, Hildeberto Jardón-Kojakhmetov
% and Sebastian Wieczorek
clear
warning off
addpath("../")

opts_fsolve = optimset('Display','off');

%% Parameters

N = 2;

ome2 = -4;

eta = 10; % adaptive parameters
alpha = pi/2; % phase shift.
kappa = 1;

Eps = 0.1;

%% Finding the equilirbrium solutions

red_Del = 100;
Del_arr = linspace(0,1,red_Del);

e_temp = [1.0307; 1.0307; 4.8576];
% e_temp = [5; 5; 4.8576];
e_mat = NaN(3,red_Del);
eig_val = NaN(3,red_Del);

for ind_Del = 1:red_Del
    Del = Del_arr(ind_Del);
    ome = [ome2 + Del; ome2];
    par = [ome; kappa; eta; alpha; Eps];

    fun = @(var) Adap_phase_osc_N(var, par, N);

    e_temp = fsolve (fun,e_temp,opts_fsolve);
    e_mat(:,ind_Del) = [mod(e_temp(1),2*pi); mod(e_temp(2),2*pi); e_temp(3)];
    % disp( fun(e_mat(:,ind_Del)) )
end
c = e_mat(1,end);

%% fix this!
opts_ode = odeset('RelTol',1e-5,'AbsTol',1e-5,'Events',...
    @(t,var)myeventfun(t,var,N));
%% simulations

phi1_init = c;
tend = 10/Eps;

res = 500; 
mu_init = 2.86;

phi1_scan = linspace(0,2*pi,res);
phi2_scan = linspace(0,2*pi,res);
tic
% define odefun
odefun = @(t, var) Adap_phase_osc_N(var, par, N);
for ind_phi1 = 1:res
    phi1_init = phi1_scan(ind_phi1);
    parfor ind_phi2 = 1:res
        initcond = [phi1_init; phi2_scan(ind_phi2); mu_init];
        [t,var] = ode45(odefun, [0 tend], initcond, opts_ode);
        if var(end,N+1) < 7
            basin_mat1(ind_phi2,ind_phi1) = 1;
        else
            basin_mat1(ind_phi2,ind_phi1) = 0;
        end
    end 
    disp(ind_phi1)
end
toc
%% simulations
% initcond = [phi1_init;5.16944 ; 0.083612];
% [t,var] = ode45(odefun, [0 1000], initcond, opts_ode);
% mu = var(:,3);
% phi1 = mod(var(:,1),2*pi);
% phi2 = mod(var(:,2),2*pi);
% figure
% plot3(mu,phi1,phi2,'.k')
% xlim([0,2*pi])
% ylim([0,2*pi])
% save('Basin_2osc.mat')
%% testing 
% initcond = [5.6; 5.6-(pi/2); -2.54545 ];
% [t,var] = ode45(odefun, [0 200], initcond,opts_test);
% figure(11)
% hold on
% plot(t, mod(var(:,1),2*pi),'k.')
% plot(t, mod(var(:,2),2*pi),'r.')
%% plotting 
% figure
% subplot(2,2,1)
% pp = pcolor(phi1_scan,phi2_scan,basin_mat1);
% pp.LineStyle = 'none';
% pp.FaceAlpha = 0.3;
% colormap ([1 1 1; 0.6 0.1 0.2])
% colormap ([1 1 1; 0 0 0; 0.5 0.5 0.5])

contour(phi1_scan, phi2_scan, basin_mat1,...
    'LineColor', [0.10,0.40,0.69],'LevelStep',1,'LineWidth',2)


%% event function
function [check,stop,direction] = myeventfun(t,var,N)
check = var(N+1)>8; %prod(var(1:N) - 40); 
stop = 1;  % Halt integration
direction = 0;
end
