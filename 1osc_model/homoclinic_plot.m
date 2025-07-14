function [gap] = homoclinic_plot(EPSS,alpha)

eta = 10; % adaptive parameters


ome = -4;

par = [ome, eta, alpha, EPSS];


% equilibria
phi_e1 = ...
    mod(asin( (ome+eta)/( sqrt( (1-eta)^2 + 2*eta*(1+cos(alpha)) ) ) ) - ...
    atan( eta*sin(alpha)/( 1+eta*cos(alpha) ) ), 2*pi) ;

mu_e1 = sin(phi_e1) - ome;

phi_e2 = mod(pi - ...
    asin( (ome+eta)/( sqrt( (1-eta)^2 + 2*eta*(1+cos(alpha)) ) ) ) - ...
    atan( eta*sin(alpha)/( 1+eta*cos(alpha) ) ), 2*pi);

mu_e2 = sin(phi_e2) - ome;



opts1 = odeset('RelTol',1e-16,'AbsTol',1e-16,'Events',@(t,var)myeventfun1(t,var));
opts2 = odeset('RelTol',1e-16,'AbsTol',1e-16,'Events',@(t,var)myeventfun2(t,var));

odefun = @(t,var)Adaptive_phase_ode(var,par);

tend =100;

% computing 

% init cond for the stable manifold of the saddel
initcond_m = [phi_e2 mu_e2] - 0.01.*[1 0];
% initcond_p = [phi_e2 mu_e2] + 0.0001.*[1 0];

[t, var_m1] = ode45(odefun,[tend 0],initcond_m,opts1);
initcond_m = [2*pi var_m1(end,2)];
[t, var_m2] = ode45(odefun,[tend 0],initcond_m,opts2);
% [t, var_p] = ode45(odefun,[tstart 0],initcond_p,opts);


phi_m1 = var_m1(:,1);
mu_m1 = var_m1(:,2);

phi_m2 = var_m2(:,1);
mu_m2 = var_m2(:,2);

figure(1);
hold on 
plot(mu_m1,phi_m1,'-k')
plot(mu_m2,phi_m2,'-k')

% ode function
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

% event function
function [check,stop,direction] = myeventfun1(t,var)
% tempvar = mod(var(1),2*pi);
check = var(1);
stop = 1;  % Halt integration
direction = 0;
end


function [check,stop,direction] = myeventfun2(t,var)
% tempvar = mod(var(1),2*pi);
check = var(1)*(var(1) - 2*pi);
stop = 1;  % Halt integration
direction = 0;
end

end 