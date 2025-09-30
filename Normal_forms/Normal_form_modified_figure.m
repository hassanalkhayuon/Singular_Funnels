% This code is created by Hassan Alkhayuon to
% simulate plot singuler basin of the modified normal form suggusted by SW
% A Research project with Serhiy Yanchuk, Hildeberto Jardón-Kojakhmetov,
% and Sebastian Wieczorek
 

clear
warning off
% Colour-blind friendly colour
[red,yellow,green,blue ] = ...
    Colour_blind_friendly_colours();
%% Parameters

a = -10;
b = 5;

EPSS = 0.1;

par = [a,b, EPSS];

%% equilibria



eq = @(xx) (tanh(a + b.*xx) + 2 - xx);

xe0 = 0;
xe1 = fzero(eq, 1);
xe2 = fzero(eq, 2);
xe3 = fzero(eq, 3);

mue0 = a;
mue1 = a + b.*xe1;
mue2 = a + b.*xe2;
mue3 = a + b.*xe3;



%% Critical manifold 

x_cm_0 =@(m)0; 
x_cm_tanh = @(m)(tanh(m)+2);


%% simulations
opts = odeset('RelTol',1e-10,'AbsTol',1e-10);

odefun = @(t,var)sing_bsin_MNF(var,par);

JJ = Jacobian_MNF([xe2;mue2],par);

[v,lambda] = eig(JJ);

if lambda(1,1)<0
    stable_dir = v(:,1);
else
    stable_dir = v(:,2);
end


initcond1 = [xe2;mue2] + 0.001.*stable_dir;
initcond2 = [xe2;mue2] - 0.001.*stable_dir;

[~,man1] = ode45(odefun, [100 0], initcond1, opts);
[~,man2] = ode45(odefun, [100 0], initcond2, opts);

x1 = man1(:,1);
mu1 = man1(:,2);
x2 = man2(:,1);
mu2 = man2(:,2);

%% potential
syms mu
potential_p = -1*int(b*tanh(mu) - mu + a + 2*b, mu);

potential_func_p = matlabFunction(potential_p);

% fplot(potential_func_p,[-10, 10],'-k','LineWidth',2)

%% shading matrix 

res_mu = 1000;
res_x = 1000;

mu_arr = linspace(-10,10,res_mu);
x_arr = linspace(0,5,res_x);
shading_mat = NaN(res_mu, res_x);

for ind_x = 1:res_x 
    x_init = x_arr(ind_x);
    parfor ind_mu = 1:res_mu
        initcond = [x_init; mu_arr(ind_mu)];
        [~,var] = ode45(odefun, [0 200], initcond, opts);
       dist = norm([var(end,2), var(end,1)] - [mue1, xe1]);

        if dist < 0.1
            shading_mat(ind_mu,ind_x) = 1;
        end
    end
    disp(ind_x)
end

%% testing 

% initcond = [0.06; 2.373];
% [~,var] = ode45(odefun, [0 200], initcond, opts);
% plot(var(:,2),var(:,1))
 
%% plotting 

figure(10); 
clf
% hold on

subplot(1,2,1)
hold on 
pp = pcolor(mu_arr, x_arr, shading_mat');

pp.LineStyle = 'none';
% 
fplot(x_cm_0,[-10 10],'--','Color',[0.7 0.7 0.7], 'LineWidth',2)
fplot(x_cm_tanh, [-10 10],'Color',[0.7 0.7 0.7], 'LineWidth',2)

% equilibria
plot(...
    mue0,xe0,'+k','MarkerSize',6,'LineWidth',2)
plot(...
    mue1,xe1,'.k','MarkerSize',20)
plot(...
    mue2,xe2,'+k','MarkerSize',6,'LineWidth',2)
plot(...
    mue3,xe3,'.k','MarkerSize',20)


plot(mu1,x1,'-r',...
    mu2,x2,'-r',...
    'LineWidth',2)

axis([-10, 10, 0, 3.5])


subplot(1,2,2)
hold on 

fplot(potential_func_p,[-10, 10],'-k','LineWidth',2)
% fplot(potential_func_m,[-3, 0],'-k', 'LineWidth',2)
% axis([-10, 10, -13, 3])

%% Typical trjectories

% init = [0.000000002, -5];
% [~,var] = ode45(odefun, [0 100], init, opts);
% x = var(:,1);
% mu = var(:,2);
% 
% figure(10); 
% hold on 
% plot(mu,x,'-k','LineWidth',1)

%% ode function
function [dvar] = sing_bsin_MNF(var,par)
%sing_bsin_MNF: the ODE function for modified normal form for singular 
% basin  
%   inputs:
%   var: system variables vector of length N+1
%   par: system parameters
%  
%   output: 
%   dvar: the time dervative for the variables


% variables
x = var(1);
mu = var(2);

% parameters
a = par(1);
b = par(2);
EPSS = par(3);

dvar(1) = x*(tanh(mu) + 2 - x);  
dvar(2) = EPSS.*(-mu + a + b.*x);
dvar = dvar';

end

%% Jacobian of the normal form

function [J] = Jacobian_MNF(var,par)
%Adap_phase_osc_N: the ODe function for adaptive phase oscillator system 
%   inputs:
%   var: system variables vector of length N+1
%   par: system parameters
%   N: number of oscillators in the network 
%   output: 
%   dvar: the time dervative for the variables


% variables
x = var(1);
mu = var(2);

% parameters
a = par(1);
b = par(2);
EPSS = par(3);

J(1,1) = tanh(mu) - 2.*x + 2; 
J(1,2) = x*(1-tanh(mu)^2); 
J(2,1) = EPSS*b;
J(2,2) = -EPSS;
end