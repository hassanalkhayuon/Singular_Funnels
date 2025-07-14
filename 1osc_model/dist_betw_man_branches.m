% This code is created by Hassan Alkhayuon to
% find the distence between the branches of stable manifold in
% terms of eps for a given rotation
% A Research project with Serhiy Yanchuk and Hildeberto Jardón-Kojakhmetov

clear
warning off
addpath('..\')

% Colour-blind friendly colour
[red,yellow,green,blue ] = ...
    Colour_blind_friendly_colours();

%% Parameters
ome = -4;
eta = 10; % adaptive parameters
kappa = 0;
alpha = pi/2; % phase shift.

rot = 1;

%% equilibria
phi_e1 = ...
    mod(asin( (ome+eta)/( sqrt( (1-eta)^2 + 2*eta*(1+cos(alpha)) ) ) ) - ...
    atan( eta*sin(alpha)/( 1+eta*cos(alpha) ) ), 2*pi) ;

mu_e1 = sin(phi_e1) - ome;

phi_e2 = mod(pi - ...
    asin( (ome+eta)/( sqrt( (1-eta)^2 + 2*eta*(1+cos(alpha)) ) ) ) - ...
    atan( eta*sin(alpha)/( 1+eta*cos(alpha) ) ), 2*pi);

mu_e2 = sin(phi_e2) - ome;


%% ode setting
res = 20;
EPS_scan = logspace(log10(0.0035),log10(0.2),res);
for ind_eps = 1:length(EPS_scan)
    EPS = EPS_scan(ind_eps);
    par = [ome, kappa, eta, alpha, EPS];
    opts = odeset(...
        'RelTol',1e-15,...
        'AbsTol',1e-15,...
        'Events',@myeventfun);
    odefun = @(t,var) Adap_phase_osc_N(var,par,1);

    TT = 100;

    % computing the stable manifold of e2
    % stable manifold of the saddel
    %
    initcond_m = [phi_e2 mu_e2] - 0.001.*[1 0];

    for ind_man = 1:rot+1

        [t, var] = ode15s(odefun,[TT 0],initcond_m,opts);

        initcond_m = [2*pi - var(end,1), var(end,2)];


        phi = var(:,1);
        mu = var(:,2);


        tstart = t(end);

    end
    end_point_1 = var(end,:);

    initcond_p = [phi_e2 mu_e2] + 0.001.*[1 0];

    for ind_man = 1:rot

        [t, var] = ode45(odefun,[TT 0],initcond_p,opts);

        initcond_p = [2*pi - var(end,1), var(end,2)];

        phi = var(:,1);
        mu = var(:,2);



        tstart = t(end);
    end

    end_point_2 = var(end,:);

    del(ind_eps) = norm(end_point_1 - end_point_2);
    
end
log_EPS = log(EPS_scan);
log_log_del = log(abs(log(del)));

%% plotting
figure(5);
hold on
cla
plot(log_EPS,log_log_del,'ok','LineWidth',2)
xlabel('$\ln(\varepsilon)$')
ylabel('$\ln(\ln(\delta))$','Rotation',0)
set(gcf, 'renderer', 'painters')
set(gca,'FontSize',15)
axis([min(log_EPS) max(log_EPS)  min(log_log_del) max(log_log_del)])
box on


%% event function
function [check,stop,direction] = myeventfun(t,var)
check = (var(1) - 2*pi)*var(1);
stop = 1;  % Halt integration
direction = 0;
end
