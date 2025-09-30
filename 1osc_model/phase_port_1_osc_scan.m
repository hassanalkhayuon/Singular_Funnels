% This code is created by Hassan Alkhayuon to
% colour the basin of attraction for adaptive phase oscillator .
% A Research project with Serhiy Yanchuk and Hildeberto Jardón-Kojakhmetov

clear
warning off

% Colour-blind friendly colour
[red,yellow,green,blue ] = ...
    Colour_blind_friendly_colours();

%% Parameters

eta = 10; % adaptive parameters
alpha = pi/2; %5.236; 1.67; 1.0472; % phase shift.

EPS = 0.01;

ome_start = -4;
% ome1 = -3.4;
% ome2 = -0.9;
% ome3 = -11;

% ome_start = ome3;

par = [ome_start, eta, alpha, EPS];


% equilibria
phi_e1 = ...
    mod(asin( (ome_start+eta)/( sqrt( (1-eta)^2 + 2*eta*(1+cos(alpha)) ) ) ) - ...
    atan( eta*sin(alpha)/( 1+eta*cos(alpha) ) ), 2*pi) ;

mu_e1 = sin(phi_e1) - ome_start;

phi_e2 = mod(pi - ...
    asin( (ome_start+eta)/( sqrt( (1-eta)^2 + 2*eta*(1+cos(alpha)) ) ) ) - ...
    atan( eta*sin(alpha)/( 1+eta*cos(alpha) ) ), 2*pi);
%%

mu_e2 = sin(phi_e2) - ome_start;

% plot(...
%     mu_e1,phi_e1,'.','MarkerSize',25,'Color',green)

%% scan

opts_noevents = odeset('RelTol',1e-10,'AbsTol',1e-10);
opts = odeset('RelTol',1e-10,'AbsTol',1e-10,'Events',@(t,var)myeventfun(t,var));
odefun = @(t,var)Adaptive_phase_ode(var,par);

%%
res_scan = 300;
phi_scan = linspace(0, 2*pi-0.01, res_scan);
mu_scan  = linspace(-4, 14, res_scan);
color_mat = zeros(res_scan,res_scan);
tol = 1e-2;

for ind_phi = 1: res_scan
    for ind_mu = 1: res_scan
        initcond_m = [phi_scan(ind_phi), mu_scan(ind_mu)];
        [t,var] = ode45(odefun,[0 50],initcond_m,opts_noevents);
        endpoint_mu = var(end,2);
        crit = abs(endpoint_mu -  mu_e1);
        if endpoint_mu <5 %crit <= tol
            color_mat(ind_phi, ind_mu) = 1;
        end
    end
    disp(ind_phi)
end




%% plotting
%

%
% basin of attraction
figure(10)
cla
hold on

pp = pcolor(mu_scan, phi_scan, color_mat);
pp.LineStyle = 'none';
pp.FaceAlpha = 0.1;
colormap ([1 1 1; 0.6 0.1 0.2])
colormap ([1 1 1; 0 0 0])
colormap([1,0,0; 1, 1, 1])
% % stable manifold of the saddel
% initcond_m = [phi_e2 mu_e2] - 0.0001.*[1 0];
% 
% rot_total = 7;
% ind_gap = 1;
% tstart = 100000;
% 
% for ind = 1:rot_total+1
% 
%     [t, var] = ode15s(odefun,[tstart 0],initcond_m,opts);
% 
%     initcond_m = [2*pi - var(end,1), var(end,2)];
%     phi = var(:,1);
%     mu = var(:,2);
% 
%     plot(mu, phi,'-','Color','k','LineWidth',2)
%     if ind == ind_gap  + 1
%         % plot(mu(end), phi(end),'.','Color',blue,'MarkerSize',20)
%     end
% end
% 
% initcond_p = [phi_e2 mu_e2] + 0.0001.*[1 0];
% 
% 
% for ind = 1:rot_total
% 
%     [t, var] = ode45(odefun,[tstart 0],initcond_p,opts);
% 
%     initcond_p = [2*pi - var(end,1), var(end,2)];
%     phi = var(:,1);
%     mu = var(:,2);
% 
%     plot(mu, phi,'-','Color','k','LineWidth',2)
%     if ind == ind_gap
%         % plot(mu(end), phi(end),'.','Color',blue,'MarkerSize',20)
%     end
% end



% critical manifold



% phi01 = linspace(0, pi/2, 1000);
% phi02 = linspace(pi/2, 3*pi/2 ,1000);
% phi03 = linspace(3*pi/2, 2*pi,1000);
% 
% mu01 = -ome_start + sin(phi01);
% mu02 = -ome_start + sin(phi02);
% mu03 = -ome_start + sin(phi03);
% 
% plot(mu01,phi01,'-','Color',[0.6 0.6 0.6])
% plot(mu02,phi02,'--','Color',[0.6 0.6 0.6])
% plot(mu03,phi03,'-','Color',[0.6 0.6 0.6])
% 
% % trajectories
% %
% 
% % cylindrical limit cycle
% %plot([mu_LC mu_LC], [0 2*pi],'-b','LineWidth',3)
% 
% % equilibria
% plot(...
%     mu_e1,phi_e1,'.k','MarkerSize',25)
% plot(...
%     mu_e2,phi_e2,'+k','MarkerSize',7,'LineWidth',2)
% 
% % plot(...
% %     mu_e1_start,phi_e1_start,'.r','MarkerSize',25)


box on
ylim([0 2*pi])
xlim([-1.2 15])
xlabel('$\mu$')
ylabel('$\varphi$')
% set(gcf, 'renderer', 'painters')
set(gca,'FontSize',15)

%% typical trajectory
% initcond = [0.0000   10.0195]; % pereodic orbit
% 
% % initcond = [6.28 9.9];
% 
% % tend =10;
% % tstart = 0;
% % figure
% % cla
% hold on
% rout_numb = 0;
% while rout_numb < 1
%     [t,var_] = ode45(odefun,[0 10],initcond,opts);
%     tend = t(end);
%     initcond = [2*pi - abs(var_(end,1)), var_(end,2)];
%     plot(var_(:,2),var_(:,1),'Color', 'k','LineWidth',2)
%     rout_numb = rout_numb + 1;
% end

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


% %% event function
function [check,stop,direction] = myeventfun(t,var)
check = (var(1) >= 2*pi) + (var(1) <= 0);
stop = 1;  % Halt integration
direction = 1;
end
