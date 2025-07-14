% This code is created by Hassan Alkhayuon to
% find basin of attraction for adaptive phase oscillator .
% A Research project with Serhiy Yanchuk and Hildeberto Jardón-Kojakhmetov

clear
warning off

% Colour-blind friendly colour
[red,yellow,green,blue ] = ...
    Colour_blind_friendly_colours();

%% Parameters

eta = 10; % adaptive parameters
alpha = -1.2; %-0.97691488076; %2.7; 2.4330401; % phase shift.

EPS = 0.04;

% addional point 
% [-1, 0.0173785;
%  -1.02, 0.0120045
%  -1.04, 0.0035001725;
% ]


ome_start = -4;
% ome1 = -3.4;
% ome2 = -0.9;
% ome3 = -11;

% ome_start = ome3;

par = [ome_start, eta, alpha, EPS];


phi_sym = sym('phi_sym');

eq_expr = -ome_start + sin(phi_sym) == eta*(1- sin(phi_sym+alpha));

phi_eq = solve(eq_expr, phi_sym,Real=true);

phi_e1 = mod(eval(phi_eq(2)),2*pi);
phi_e2 = mod(eval(phi_eq(1)),2*pi);

mu_e1 = sin(phi_e1) - ome_start;
mu_e2 = sin(phi_e2) - ome_start;

%% ode functions and options

opts_noevents = odeset('RelTol',1e-10,'AbsTol',1e-10);
opts = odeset('RelTol',1e-10,'AbsTol',1e-10,'Events',@(t,var)myeventfun(t,var));

odefun = @(t,var)Adaptive_phase_ode(var,par);
% Jacfun = @(var)Adaptive_phase_ode(var,par);

% e1 = [phi_e1, mu_e1];
% e2 = [phi_e2, mu_e2];
% 
% jj = MyJacobian(Jacfun,e1);
% 
% lambda = eig(jj);
% 
% re_lambda = real(lambda(1));
% 
% figure(6)
% plot(alpha,re_lambda)


%% plotting
%


% basin of attraction

figure(4);
cla
hold on

% critical manifold

phi01 = linspace(0, pi/2, 1000);
phi02 = linspace(pi/2, 3*pi/2 ,1000);
phi03 = linspace(3*pi/2, 2*pi,1000);

mu01 = -ome_start + sin(phi01);
mu02 = -ome_start + sin(phi02);
mu03 = -ome_start + sin(phi03);

phi_sap = linspace(0,2*pi,100);
mu_sap = eta*(1 - sin(phi_sap + alpha)); 

plot(mu01,phi01,'-','Color',[0.6 0.6 0.6])
plot(mu02,phi02,'--','Color',[0.6 0.6 0.6])
plot(mu03,phi03,'-','Color',[0.6 0.6 0.6])

% stable manifold of the saddel
initcond_m = [phi_e2 mu_e2] - 0.001.*[0 1];

rot_total = 200;
tstart = 100;

for ind = 1:rot_total+1

    [t, var] = ode15s(odefun,[tstart 0],initcond_m,opts);

    initcond_m = [mod(var(end,1),2*pi), var(end,2)];
    phi = var(:,1);
    mu = var(:,2);
    plot(mu, phi,'-','Color',[0.00,0.45,0.74],'LineWidth',2)
end

initcond_p = [phi_e2 mu_e2] + 0.0001.*[0 1];


for ind = 1:rot_total

    [t, var] = ode15s(odefun,[tstart 0],initcond_p,opts);
    
    initcond_p = [mod(var(end,1),2*pi), var(end,2)];
    phi = var(:,1);
    mu = var(:,2);

    plot(mu, phi,'-','Color',[0.00,0.45,0.74],'LineWidth',2)
end

% unstable manifold of the saddel
initcond_m = [phi_e2 mu_e2] - 0.001.*[0 1];

rot_total = 1000;
tstart = 100;

for ind = 1:rot_total+1

    [t, var] = ode15s(odefun,[0 tstart],initcond_m,opts);

    initcond_m = [mod(var(end,1),2*pi), var(end,2)];
    phi = var(:,1);
    mu = var(:,2);
    plot(mu, phi,'-k','LineWidth',1)
end

initcond_p = [phi_e2 mu_e2] + 0.0001.*[0 1];


for ind = 1:rot_total

    [t, var] = ode15s(odefun,[0 tstart],initcond_p,opts);

    initcond_p = [mod(var(end,1),2*pi), var(end,2)];
    phi = var(:,1);
    mu = var(:,2);

    plot(mu, phi,'-k','LineWidth',1)
end


% trajectories
%

% cylindrical limit cycle

[t, var_LC] = ode15s(odefun,[0 10],[0 10],opts);
plot(var_LC(:,2),var_LC(:,1),'-k','LineWidth',3)

% unstable cylindrical limit cycle
% [t, var_LC] = ode15s(odefun,[10 0],[2*pi 5.7],opts);
% plot(var_LC(:,2),var_LC(:,1),'-k','LineWidth',3)


% equilibria
plot(...
    mu_e1,phi_e1,'.w','MarkerSize',25)
plot(...
    mu_e1,phi_e1,'ok','MarkerSize',6,'LineWidth',2)
plot(...
    mu_e2,phi_e2,'.w','MarkerSize',25)
plot(...
    mu_e2,phi_e2,'ok','MarkerSize',6,'LineWidth',2)


box on
% axis([0 12 0 2*pi])
% xlabel('$\mu$')
% ylabel('$\varphi$')
% % set(gcf, 'renderer', 'painters')
% set(gca,'FontSize',15)


%% additional trajectories trajectory
% initcond = [phi_e1 - 0.001, mu_e1];
% 
% 
% t0 = 0;
% t1 = 1000;
% 
% figure(4);
% hold on
% while t0 < t1
%     [t,var_] = ode15s(odefun,[t0 t1],initcond,opts);
%     t0 = t(end);
%     initcond = [mod(var_(end,1),2*pi), var_(end,2)];
%     plot(var_(:,2),var_(:,1),'Color', 'b')
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


%% event function
function [check,stop,direction] = myeventfun(t,var)
check = (var(1) >= 2*pi) + (var(1) <= 0);
stop = 1;  % Halt integration
direction = 1;
end

%% numarical jacobian (From Chat GPT, need to be checked )

function J = MyJacobian(func, x, h)
    % computeJacobianNumeric computes the Jacobian matrix of a given function
    % using numerical differentiation.
    %
    % Inputs:
    %   func - a function handle representing the vector-valued function
    %   x    - a vector of input values at which to evaluate the Jacobian
    %   h    - a small step size for finite differences (optional, default is 1e-6)
    %
    % Output:
    %   J - the Jacobian matrix

    % Set default step size if not provided
    if nargin < 3
        h = 1e-6; % Default step size
    end

    % Number of input variables
    n = length(x);
    
    % Evaluate the function at the original point
    f0 = func(x);
    
    % Number of output functions
    m = length(f0);
    
    % Initialize the Jacobian matrix
    J = zeros(m, n);
    
    % Compute the Jacobian using finite differences
    for i = 1:n
        % Create a perturbation vector
        perturb = zeros(n, 1);
        perturb(i) = h;
        
        % Evaluate the function at x + h and x - h
        f_plus = func(x + perturb);
        f_minus = func(x - perturb);
        
        % Compute the finite difference
        J(:, i) = (f_plus - f_minus) / (2 * h);
    end
end