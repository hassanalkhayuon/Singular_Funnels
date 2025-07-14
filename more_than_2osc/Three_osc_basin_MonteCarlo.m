% This code is created by Hassan Alkhayuon to
% compute the size of basin of attraction for 3 coupled adaptive phase
% oscillators by monte carlo
% A Research project with Serhiy Yanchuk,  Hildeberto Jardón-Kojakhmetov,
% and Sebastian Wieczorek

clear
warning off


%% Parameters

N = 3;

ome = linspace(-3,-4,N)';


eta = 10; % adaptive parameters
alpha = pi/2; % phase shift.
kappa = 1;


opts = odeset('RelTol',1e-5,'AbsTol',1e-5,...
    'Events',@(t,var)myeventfun(t,var,N));


%% simulations

EPS_res = 25;
EPS_arr = linspace(1e-2,1e-1,EPS_res);

% the volume of the periodic orbit 
basin_volume = zeros(size(EPS_arr));


number_of_simulations = 100000;
% rng(0)

% figure(1);
% hold on
% cla
for ind_EPS = 1:EPS_res
    EPSS = EPS_arr(ind_EPS);
    tend = 10/EPSS;
    basin_vol = 0;

    par = [ome;  kappa; eta; alpha; EPSS];
    odefun = @(t,var)Adap_phase_osc_N(var,par,N);

    parfor ind = 1:number_of_simulations

        rr = rand(1,4);

        phi1_init = rr(1).*2.*pi; 
        phi2_init = rr(2).*2.*pi;
        phi3_init = rr(3).*2.*pi;
        mu_init = rr(4)*4 - 4;

        initcond = [phi1_init; phi2_init; phi3_init; mu_init];

        % initcond = [4.1033; phi2_init; 1];

        
        [t,var] = ode45(odefun, [0 tend], initcond, opts);

        % endpoint  = [mod(var(end,1:N),2*pi), var(end,N+1)];
        % figure(13);
        % hold on
        % plot(t,var(:,N+1))
        % drawnow 
        if var(end,N+1) > 7
            basin_vol = basin_vol + 1;
        end

    end
    basin_volume(ind_EPS) = basin_vol/number_of_simulations;
    disp(ind_EPS)
    
end
save("MC_basin_vul_2osc_unsync_1mil.mat")

%% plotting
figure(14);
hold on 

plot(1./EPS_arr,log(basin_volume),'.-','MarkerSize',20)


%% event function
% function [check,stop,direction] = myeventfun(t,var,N)
% check = var(N+1)>8; %prod(var(1:N) - 40);
% stop = 1;  % Halt integration
% direction = 0;
% end

function [check, stop, direction] = myeventfun(t, var, N)
    persistent lastValue; % Persistent variable to store the last value
    threshold = 1e-3;    % Define a threshold for change

    if isempty(lastValue)
        lastValue = var(N+1); % Initialize on the first call
    end

    % Check if the change is less than the threshold
    change = norm(var(N+1) - lastValue);
    check = change < threshold; % Event occurs when change is small

    % Update lastValue for the next call
    lastValue = var(N+1);

    stop = 1;  % Halt integration
    direction = 0; % No specific direction for the event
end