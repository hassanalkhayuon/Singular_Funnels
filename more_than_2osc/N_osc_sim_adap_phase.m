% This code is created by Hassan Alkhayuon to
% simulate adaptive phase oscillator.
% A Research project with Serhiy Yanchuk

clear
warning off
eps = 1e-10;
%% Parameters

N = 2;

ome = [-4; -4];
initcond = [1; 5; 5];

eta = 10; % adaptive parameters
alpha = pi/2; % phase shift.
kappa = 0.1;

EPSS = 0.1;

par = [ome;  kappa; eta; alpha; EPSS];

%% equilibria only for 1 d
phi_e1 = ...
    mod( asin( (ome+eta)/( sqrt( (1-eta)^2 + 2*eta*(1+cos(alpha)) ) ) ) - ...
    atan( eta*sin(alpha)/( 1+eta*cos(alpha) ) ), 2*pi) ;

mu_e1 = sin(phi_e1) - ome;

phi_e2 = mod(pi - ...
    asin( (ome+eta)/( sqrt( (1-eta)^2 + 2*eta*(1+cos(alpha)) ) ) ) - ...
    atan( eta*sin(alpha)/( 1+eta*cos(alpha) ) ), 2*pi);

mu_e2 = sin(phi_e2) - ome;


%% simulations

opts = odeset('RelTol',1e-10,'AbsTol',1e-10,'Events',@(t,var)myeventfun(t,var,N));
odefun = @(t,var)Adap_phase_osc_N(var,par,N);
tstart = 0;
tend = 1000;
figure(13)
hold on 
plot3(initcond(3),initcond(1),initcond(2),'.k','MarkerSize',20)
while tstart < tend
    [t,var] = ode45(odefun, [tstart tend], initcond, opts);
    tstart = t(end);
    phi = var(:,1:N);
    mu  = var(:,N+1);
    initcond = [(mod(phi(end,:),2*pi))'; mu(end)];
    plot3(mu,phi(:,1),phi(:,2),'k')
end
box on
xlabel('$\mu$')
ylabel('$\varphi_1$')
zlabel('$\varphi_2$','Rotation',0)
xlim([0 12])
ylim([0 2*pi])
zlim([0 2*pi])
%%
% tend = 100;
% tstart = 0;
% figure(10);
% hold on 
% while tstart < tend
%     [t,var] = ode45(odefun,[tstart tend],initcond,opts);
%     tstart = t(end);
%     initcond = var(end,:);
%     plot(var(:,2),var(:,1),'')
% end

% mu_LC = sol(end).y(2,end);





%% event function
function [check,stop,direction] = myeventfun(t,var,N)
check = prod( var(1:N)-2*pi).*prod(var(1:N));
stop = 1;  % Halt integration
direction = 0;
end
