% This code is created by Hassan Alkhayuon to
% simulate plot singuler basin of the normal form 
% A Research project with Serhiy Yanchuk and Hildeberto Jardón-Kojakhmetov
 

clear
warning off
% Colour-blind friendly colour
addpath('../')
[red,yellow,green,blue ] = ...
    Colour_blind_friendly_colours();
%% Parameters

a = -1;
b = 2.3;

EPSS = 0.1;

par = [a,b,EPSS];

%% equilibria

me_p = ( (b.^2 + 2.*a) + sqrt( (2.*a+b.^2)^2 - 4*a^2) ) ./ 2;
me_m = ( (b.^2 + 2.*a) - sqrt( (2.*a+b.^2)^2 - 4*a^2) ) ./ 2;
e1 = [a; 0];
e2 = [me_p; sqrt(me_p)];
e3 = [me_m; sqrt(me_m)];

%% Critical manifold 

x_cm_0 =@(m)0; 
x_cm_sqrt = @(m)sqrt(m)*(m>=0);



%% simulations
opts = odeset('RelTol',1e-10,'AbsTol',1e-10);

odefun = @(t,var)sing_bsin_NF(var,par);

JJ = Jacobian_NF([e3(2);e3(1)],par);

[v,lambda] = eig(JJ);

if lambda(1,1)<0
    stable_dir = v(:,1);
else
    stable_dir = v(:,2);
end


initcond1 = [e3(2);e3(1)] + 0.001.*stable_dir;
initcond2 = [e3(2);e3(1)] - 0.001.*stable_dir;

[~,man1] = ode45(odefun, [100 0], initcond1, opts);
[~,man2] = ode45(odefun, [100 0], initcond2, opts);

x1 = man1(:,1);
mu1 = man1(:,2);
x2 = man2(:,1);
mu2 = man2(:,2);

%% potential
syms mu 
potential_p = int(mu - a - b*sqrt(mu), mu);

potential_func_p = matlabFunction(potential_p);

potential_func_m = @(mu)( (mu^2/2) - a*mu);

%% shading matrix 

res_mu = 10;
res_x = 10;

mu_arr = linspace(-5,5,res_mu);
x_arr = linspace(-5,5,res_x);
shading_mat = NaN(res_mu, res_x);

for ind_x = 1:res_x 
    x_init = x_arr(ind_x);
    parfor ind_mu = 1:res_mu
        initcond = [x_init; mu_arr(ind_mu)];
        [~,var] = ode45(odefun, [0 200], initcond, opts);
       dist = norm([var(end,2), var(end,1)] - e2');

        if dist > 0.1
            shading_mat(ind_mu,ind_x) = 1;
        end
    end
    disp(ind_x)
end

%% testing 

initcond = [0.06; 2.373];
[~,var] = ode45(odefun, [0 200], initcond, opts);
plot(var(:,2),var(:,1))
 
%% plotting 

figure(10); 
clf
subplot(1,2,1)
hold on 
pp = pcolor(mu_arr, x_arr, shading_mat');

pp.LineStyle = 'none';

fplot(x_cm_0,[0 10],'--','Color',[0.7 0.7 0.7], 'LineWidth',2)
fplot(x_cm_sqrt, [-10 10],'Color',[0.7 0.7 0.7], 'LineWidth',2)


plot(mu1,x1,'-r',...
    mu2,x2,'-r',...
    'LineWidth',2)
axis([-5, 5, 0, 2.5])

plot(...
    e1(1),e1(2),'.k',...
    e2(1),e2(2),'.k',...
    'MarkerSize', 20)

plot(...
    e3(1),e3(2),'+k',...
    'LineWidth',2, 'MarkerSize',6)


subplot(1,2,2)
hold on 

fplot(potential_func_p,[0, 7],'-k','LineWidth',2)
fplot(potential_func_m,[-3, 0],'-k', 'LineWidth',2)
axis([-2.3, 5, -0.7, .4])

%% Typical trjectories

init = [0.02, 3.9];
[~,var] = ode45(odefun, [0 300], init, opts);
x = var(:,1);
mu = var(:,2);

figure(1); plot(mu,x,'-k','LineWidth',1)

%% ode function
function [dvar] = sing_bsin_NF(var,par)
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

dvar(1) = x*(mu-x.^2);  
dvar(2) = EPSS.*(-mu + a + b.*x);
dvar = dvar';

end

%% Jacobian of the normal form

function [J] = Jacobian_NF(var,par)
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

J(1,1) = mu - 3*x^2; 
J(1,2) = x; 
J(2,1) = EPSS*b;
J(2,2) = -EPSS;
end