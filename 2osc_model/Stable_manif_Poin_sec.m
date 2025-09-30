% This code is created by Hassan Alkhayuon to
% a parameter plane for the 2 coupled adaptive phase
% oscillators to show the stable manifold of the saddel poin on
% Poincaré Section
% A Research project with Serhiy Yanchuk, Hildeberto Jardón-Kojakhmetov
% and Sebastian Wieczorek

clear
clc
warning off
addpath("../")

options = optimset('Display','off');
%% Parameters
N = 2;
eta = 10; % adaptive parameters
alpha = pi/2; % phase shift
kappa = 1;
ome2 = -4;
Eps = 0.05;

%% Finding the equilirbrium solutions

res_Del = 100;
Del_arr = linspace(0,1,res_Del);

e_temp = [1.0307; 1.0307; 4.8576];
e_temp = [5; 5; 4.8576];
e_mat = NaN(3,res_Del);
eig_val = NaN(3,res_Del);

for ind_Del = 1:res_Del
    Del = Del_arr(ind_Del);
    ome = [ome2 + Del; ome2];
    par = [ome; kappa; eta; alpha; Eps];

    fun = @(var) Adap_phase_osc_N(var, par, N);

    e_temp = fsolve (fun,e_temp,options);
    e_mat(:,ind_Del) = [mod(e_temp(1),2*pi); mod(e_temp(2),2*pi); e_temp(3)];
    % disp( fun(e_mat(:,ind_Del)) )
end
%% informations about the Poincaré Section
c = e_mat(1,end);
opts = odeset('RelTol',1e-5,'AbsTol',1e-5,'Events',...
    @(t,var)event_fun(t, var, c));

%% unstable manifold
JJ = MyJacobian(fun,e_mat(:,end));
[eig_vic, eig_val] = eig(JJ);



v1 = eig_vic(:,1);
v2 = eig_vic(:,2);
v3 = eig_vic(:,3);

res_man = 1000;
figure(13);
hold on
for ind_man = 1:res_man
    a = rand - 0.5;
    b = rand - 0.5;
    initconds = e_mat(:,end) + 1e-3*a*v2 + 1e-5*b*v3;

    odefun = @(t,var) Adap_phase_osc_N(var, par, N);
    [t,var,te, vare, ie] = ode45(odefun, [100 0], initconds,opts);
    mu = vare(:,3);
    phi1 = mod(vare(:,1),2*pi);
    phi2 = mod(vare(:,2),2*pi);
    if abs(phi1-c)<=0.001
        plot(mu, phi2,'.b')
        drawnow
        disp(ind_man)
    end
end



plot(e_mat(3,end), e_mat(2,end),'.r','MarkerSize',30)
%% Numerical Jacobian
function J = MyJacobian(fun, x, epsilon)
if nargin < 3
    epsilon = 1e-8; % Default value for epsilon
end

% Number of variables
n = length(x);
% Number of functions
m = length(fun(x));

% Initialize Jacobian matrix
J = zeros(m, n);

% Compute the Jacobian using finite differences
for ind = 1:n
    % Create a perturbed vector
    x_p_eps = x;
    x_p_eps(ind) = x_p_eps(ind) + epsilon;

    % Evaluate the function at the perturbed point
    J(:, ind) = (fun(x_p_eps) - fun(x)) / epsilon;
end
end
%%
function [value, isterminal, direction] = event_fun(t, var, c)
phi1 = mod(var(1),2*pi);
phi2 = mod(var(2),2*pi);
value = phi1 - c; 
isterminal = 0;  % don't stop the integration 
direction = 0;   % direction of crossing
end