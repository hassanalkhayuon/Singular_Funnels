% This code is created by Hassan Alkhayuon to
% compute basin of attraction for 2 coupled adaptive phase oscillators.
% A Research project with Serhiy Yanchuk and Hildeberto Jardón-Kojakhmetov

clear
warning off


%% Parameters

N = 2;

ome = [-4; -3];


eta = 10; % adaptive parameters
alpha = pi/2; % phase shift.
kappa = 1;

EPSS = 0.05;

par = [ome;  kappa; eta; alpha; EPSS];

odefun = @(t,var)Adap_phase_osc_N(var,par,N);

opts = odeset('RelTol',1e-10,'AbsTol',1e-10,...
    'Events',@(t,var)myeventfun(t,var,N));
opts_test = odeset('RelTol',1e-10,'AbsTol',1e-10);


%% simulations

phi2_init = 0.068793860248536;
tend = 1000;

res = 100; 
mu_scan = linspace(-4,0,res);
phi1_scan = linspace(0,2*pi,res);
tic
for ind_mu = 1:res
    mu_init = mu_scan(ind_mu);
    parfor ind_phi1 = 1:res
        initcond = [phi1_scan(ind_phi1); phi2_init; mu_init];
        [t,var] = ode45(odefun, [0 tend], initcond, opts);
        if var(end,N+1) < 5
            basin_mat1(ind_phi1,ind_mu) = 1;
        else
            basin_mat1(ind_phi1,ind_mu) = 0;
        end
    end 
    disp(ind_mu)
end
toc
% save('Basin_2osc.mat')
%% testing 
initcond = [5.6; 5.6-(pi/2); -2.54545 ];
[t,var] = ode45(odefun, [0 200], initcond,opts_test);
figure(11)
hold on
plot(t, mod(var(:,1),2*pi),'k.')
plot(t, mod(var(:,2),2*pi),'r.')
%% plotting 
figure(15)
% subplot(2,2,1)
pp = pcolor(mu_scan,phi1_scan,basin_mat1);
pp.LineStyle = 'none';
pp.FaceAlpha = 1;
% colormap ([1 1 1; 0.6 0.1 0.2])
colormap ([1 1 1; 0 0 0; 0.5 0.5 0.5])


%% event function
function [check,stop,direction] = myeventfun(t,var,N)
check = var(N+1)>6; %prod(var(1:N) - 40); 
stop = 1;  % Halt integration
direction = 0;
end
