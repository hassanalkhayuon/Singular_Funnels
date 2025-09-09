% This code is created by Hassan Alkhayuon to
% compute basin of attraction for 2 coupled adaptive phase oscillators.
% A Research project with Serhiy Yanchuk and Hildeberto Jardón-Kojakhmetov

clear
warning off
addpath("../")

opts_fsolve = optimset('Display','off');

%% Parameters

N = 2;

ome2 = -4;

eta = 10; % adaptive parameters
alpha = pi/2; % phase shift.
kappa = 1;

Eps = 0.07;

%% Finding the equilirbrium solutions

red_Del = 100;
Del_arr = linspace(0,1.5,red_Del);

e_temp = [1.0307; 1.0307; 4.8576];
e_mat = NaN(3,red_Del);
eig_val = NaN(3,red_Del);

for ind_Del = 1:red_Del
    Del = Del_arr(ind_Del);
    ome = [ome2 + Del; ome2];
    par = [ome; kappa; eta; alpha; Eps];

    fun = @(var) Adap_phase_osc_N(var, par, N);

    e_temp = fsolve (fun,e_temp,opts_fsolve);
    e_mat(:,ind_Del) = [mod(e_temp(1),2*pi); mod(e_temp(2),2*pi); e_temp(3)];
    % disp( fun(e_mat(:,ind_Del)) )
end
c = e_mat(1,end);

%% fix this!
opts_ode = odeset('RelTol',1e-5,'AbsTol',1e-5,'Events',...
    @(t,var)myeventfun(t,var,N));
%% simulations

tend = 10/Eps;

res = 150; 

mu_scan = linspace(4.7,-1,100); %e_mat(3,end);



phi1_scan = linspace(0,2*pi,res);
phi2_scan = linspace(0,2*pi,res);
tic
% define odefun
odefun = @(t, var) Adap_phase_osc_N(var, par, N);

V = VideoWriter('basin_movie_del_1p5.avi'); % choose filename
V.FrameRate = 5; % frames per second
open(V);

for ind_m = 1:100
    mu_init = mu_scan(ind_m);
    for ind_1 = 1:res
        phi1_in = phi1_scan(ind_1);
        parfor ind_2 = 1:res
            initcond = [phi1_in; phi2_scan(ind_2); mu_init];
            [t,var] = ode45(odefun, [0 tend], initcond, opts_ode);
            if var(end,N+1) < 7
                basin_mat1(ind_2,ind_1) = 1;
            else
                basin_mat1(ind_2,ind_1) = 0;
            end
        end
    end
    disp(ind_m)
    figure(1)
    cla
    pp = pcolor(phi1_scan,phi2_scan,basin_mat1);
    pp.LineStyle = 'none';
    pp.FaceAlpha = 1;
    colormap ([1 1 1; 0 0 0; 0.5 0.5 0.5])
    xlabel('$\varphi_1$')
    ylabel('$\varphi_2$','Rotation',0)
    title(['$\mu = $', num2str(mu_scan(ind_m))])
    clim([0 1])
    drawnow

    frame = getframe(gcf);
    writeVideo(V, frame);
end
close(V)
toc

%% event function
function [check,stop,direction] = myeventfun(t,var,N)
check = var(N+1)>8; %prod(var(1:N) - 40); 
stop = 1;  % Halt integration
direction = 0;
end
