% This code is created by Hassan Alkhayuon to
% compute the size of basin of attraction for 2 coupled adaptive phase
% oscillators by monte carlo on varying Delta_mu
% A Research project with Serhiy Yanchuk and Hildeberto Jardón-Kojakhmetov

clear
warning off


%% Parameters

N = 2;

eta = 10; % adaptive parameters
alpha = pi/2; % phase shift.
kappa = 1;
ome2 = -4;

opts = odeset('RelTol',1e-5,'AbsTol',1e-5,...
    'Events',@(t,var)myeventfun(t,var,N));


%% simulations

Eps = 0.07;
tend = 150;

res_Del = 100;

Del_arr = linspace(0,5,res_Del);
basin_volume_l = zeros(size(Del_arr));
basin_volume_e = zeros(size(Del_arr));
basin_volume_n = zeros(size(Del_arr));

number_of_simulations = 100000;

rr = rand(number_of_simulations, 3);

tic
for ind_Del = 1:res_Del
    
    basin_vol_l = 0;
    basin_vol_e = 0;
    basin_vol_n = 0;
    
    ome1 = ome2 + Del_arr(ind_Del);
    ome = [ome1;ome2];
    par = [ome;  kappa; eta; alpha; Eps];
    odefun = @(t,var)Adap_phase_osc_N(var,par,N);

    parfor ind = 1:number_of_simulations

        phi1_init = rr(ind, 1) * 2 * pi; 
        phi2_init = rr(ind, 2) * 2 * pi;
        mu_init = rr(ind, 3) * 10 - 10;

        initcond = [phi1_init; phi2_init; mu_init];

        % initcond = [4.1033; phi2_init; 1];

        
        [t,var] = ode15s(odefun, [0 tend], initcond, opts);

        endpoint  = [mod(var(end,1:N),2*pi), var(end,N+1)];

        if var(end,N+1) > 7
            basin_vol_l = basin_vol_l + 1;
        elseif abs(var(end,N+1) - 2.864) <= 1e-3
            basin_vol_e = basin_vol_e + 1;
        else
            basin_vol_n = basin_vol_n + 1;
        end
    end
    basin_volume_l(ind_Del) = basin_vol_l/number_of_simulations;
    basin_volume_e(ind_Del) = basin_vol_e/number_of_simulations;
    basin_volume_n(ind_Del) = basin_vol_n/number_of_simulations;

    disp(res_Del - ind_Del)
    
end
save("Basin_vul_2osc_eps0p06_vary_ome.mat")
toc
%% plotting
figure(8);
clf
plot(Del_arr,basin_volume_l,'.-m','MarkerSize',20)


%% event function

function [check, stop, direction] = myeventfun(t, var, N)
    check = var(N+1)>8; 
    stop = 1;
    direction = 0; 
end