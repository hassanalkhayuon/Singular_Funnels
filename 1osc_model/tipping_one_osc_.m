% This code is created by Hassan Alkhayuon to
% simulate tipping in one adaptive phase oscillator 

clear
warning off
% Colour-blind friendly colour
[red,yellow,green,blue ] = ...
    Colour_blind_friendly_colours();
%% Parameters

N = 1;

color_plot = green; [0.65 0.65 0.65];

% ome_end = -6.5;
ome_start = -3.4; ome_del = -3.1;
ome_start = -0.9; ome_del = -5.6;
ome_start = -11; ome_del = 4.5;

% ome_start = -0.9;
% ome_del = -5.6;
rate = 10;   % fast (tipping) 10, slow (no tipping) 0.1

eta = 10; % adaptive parameters
alpha = pi/2; % phase shift.

EPSS = 0.2;

par = [ome_start; ome_del; rate; eta; alpha; EPSS];

%% equilibria only for 1 d
phi_e1 = ...
    mod( asin( (ome_start+eta)/( sqrt( (1-eta)^2 + 2*eta*(1+cos(alpha)) ) ) ) - ...
    atan( eta*sin(alpha)/( 1+eta*cos(alpha) ) ), 2*pi) ;

mu_e1 = sin(phi_e1) - ome_start;

phi_e2 = mod(pi - ...
    asin( (ome_start+eta)/( sqrt( (1-eta)^2 + 2*eta*(1+cos(alpha)) ) ) ) - ...
    atan( eta*sin(alpha)/( 1+eta*cos(alpha) ) ), 2*pi);

mu_e2 = sin(phi_e2) - ome_start;

initcond = [phi_e1; mu_e1];

%% simulations

opts = odeset('RelTol',1e-10,'AbsTol',1e-10);
odefun = @(t,var)Adap_phase_osc_Nonaut(t,var,par);
tstart = -1000;
tend = 1000;

[t,var] = ode45(odefun, [tstart tend], initcond, opts);
tstart = t(end);
phi = var(:,1);
mu  = var(:,2);
ome_t = ome_start + (ome_del./2).*(tanh(rate.*t) + 1);

%% plotting 

figure(10); 

subplot(2,1,1)
hold on 
% plot(t,phi,'Color',green)
box on
set(gcf, 'renderer', 'painters')
set(gca,'FontSize',15)

subplot(2,1,2)
hold on 
plot(t,ome_t,'Color',green)
box on
set(gcf, 'renderer', 'painters')
set(gca,'FontSize',15)
xlabel('$t$')
%% phase figure 
opts = odeset('RelTol',1e-10,'AbsTol',1e-10,'Events',@(t,var)myeventfun(t,var));
tend =100;
tstart = -100;

initcond = [phi_e1; mu_e1];

% figure(3);
% cla
hold on
while tstart < tend
    [t,var_] = ode45(odefun,[tstart tend],initcond,opts);
    phi = var_(:,1);
    mu  = var_(:,2);
    tstart = t(end);
    initcond = [2*pi - abs(var_(end,1)), var_(end,2)];
    % subplot(2,1,1)
    hold on
    plot(t,phi,'Color',color_plot,'LineWidth',1.5)
    box on
    % set(gcf, 'renderer', 'painters')
    set(gca,'FontSize',15)
end

%% event function
function [check,stop,direction] = myeventfun(t,var)
check = (var(1) >= 2*pi) + (var(1) <= 0);
stop = 1;  % Halt integration
direction = 1;
end
