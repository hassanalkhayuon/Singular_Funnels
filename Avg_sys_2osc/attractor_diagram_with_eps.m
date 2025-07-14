% This code is created by Hassan Alkhayuon to compute an attractor diagram
% for the 2 phase oscelator network
% This is a Research project with Serhiy Yanchuk,
% Hildeberto Jardón-Kojakhmetov, and Sebastian Wieczorek
clear
warning off
addpath('../')

%%  parameters
N = 2;

ome = [-3; -4];


eta = 10; % adaptive parameters
alpha = pi/2; % phase shift.
kappa = 1;

EPS_arr = linspace(0.01,0.1,100);

%% ode set up
opts = odeset('RelTol',1e-5,'AbsTol',1e-5);

tspan1 = [0 1000];
tspan2 = [0 100];

% init conditions
phi1_init = 5.9;
phi2_init = 5.24;
mu_init = 2.8;

initcond = [phi1_init; phi2_init; mu_init];
    %%
figure(11);
clf
hold on
axis([EPS_arr(1), EPS_arr(end),0,12])

% figure(12);
% clf
% hold on
% axis([EPS_arr(1), EPS_arr(end),0,2*pi])
% 
% figure(13);
% clf
% hold on
% axis([EPS_arr(1), EPS_arr(end),-3, 15])

for ind_eps = 1:length(EPS_arr)
    EPS = EPS_arr(ind_eps);
    par = [ome;  kappa; eta; alpha; EPS];
    odefun = @(t,var)Adap_phase_osc_N(var,par,N);
    [~,var_temp] = ode15s(odefun, tspan1, initcond, opts);
    new_init = var_temp(end,:);
    [t,var] = ode15s(odefun, tspan2, new_init, opts);

    phi1 = mod(var(:,1),2*pi);
    phi2 = mod(var(:,2),2*pi);
    mu = var(:,end);

    % phi1 local max
    phi1_max = phi1(islocalmax(phi1));
    if isempty(phi1_max) 
        phi1_max = phi1(end);
    end

    % phi1 local min
    phi1_min = phi1(islocalmin(phi1));
    if isempty(phi1_min)
        phi1_min = phi1(end);
    end

    % phi2 local max 
    phi2_max = phi2(islocalmax(phi2));
    if isempty(phi2_max)
        phi2_max = phi2(end);
    end

    % phi2 local min
    phi2_min = phi2(islocalmin(phi2));
    if isempty(phi2_min)
        phi2_min = phi2(end);
    end

    % mu max
    mu_max = mu(islocalmax(mu));
    if isempty(mu_max)
        mu_max = mu(end);
    end

    % mu min
    mu_min = mu(islocalmin(mu));
    if isempty(mu_min)
        mu_min = mu(end);
    end

    figure(11);
    plot(EPS,phi1_max,'.r')
    plot(EPS,phi1_min,'.r')

    % figure(12);
    plot(EPS,phi2_max,'.b')
    plot(EPS,phi2_min,'.b')
    % 
    % figure(13);
    plot(EPS,mu_max,'.k',...
        EPS,mu_min,'.k')


    drawnow()
end
% plot(time,phi1,'-k',...
%     time(ind_min),phi1(ind_min),'*r',...
%     time(ind_max),phi1(ind_max),'*r')
%% additionla plotting 

EPS = 0.06925;
par = [ome;  kappa; eta; alpha; EPS];
odefun = @(t,var)Adap_phase_osc_N(var,par,N);

% [~,var_temp] = ode15s(odefun, tspan1, initcond, opts);
% new_init = var_temp(end,:);
[t,var] = ode15s(odefun, [0 10], initcond, opts);

phi1 = mod(var(:,1),2*pi);
phi2 = mod(var(:,2),2*pi);
mu = var(:,end);
figure(15)
hold
plot(mu,phi1)