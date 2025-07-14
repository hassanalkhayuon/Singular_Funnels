% This code is created by Hassan Alkhayuon to
% compute the size of basin of attraction for 1 coupled adaptive phase
% oscillators by monte carlo
% A Research project with Serhiy Yanchuk and Hildeberto Jardón-Kojakhmetov

clear
warning off
addpath('..\')

%% Parameters

N = 1;

ome = -4;


eta = 10; % adaptive parameters
alpha = pi/2; % phase shift.
kappa = 0;

opts = odeset('RelTol',1e-10,'AbsTol',1e-10,...
    'Events',@(t,var)myeventfun(t,var,N));


%% simulations

EPS_res = 50;
EPS_arr = linspace(1e-3,1e-1,EPS_res);

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

        rr = rand(1,3);

        phi1_init = rr(1).*2.*pi; 
        phi2_init = rr(2).*2.*pi; 
        mu_init = rr(3)*10 - 200; 

        initcond = [phi1_init; mu_init];

        % initcond = [4.1033; phi2_init; 1];

        
        [t,var] = ode45(odefun, [0 tend], initcond, opts);

        endpoint  = [mod(var(end,1:N),2*pi), var(end,N+1)];

        if var(end,N+1) > 4
            basin_vol = basin_vol + 1;
        %     plot(t,var(:,end),'r')
        %     xlim([0 50])
        %     ylim([0 11])
        %     drawnow()
        % else
        %     plot(t,var(:,end),'k')
        %     xlim([0 50])
        %     ylim([0 11])
        %     drawnow()
        end

    end
    basin_volume(ind_EPS) = basin_vol/number_of_simulations;
    disp(ind_EPS)
end


%% plotting

figure(20);
hold on 

valume_index_cut = find(basin_volume>1e-4,1);
delta = 1./EPS_arr(valume_index_cut:end);

yy = log(basin_volume(valume_index_cut:end)) - log(delta);

P1 = polyfit(delta,yy,1);

B = P1(2);
A = P1(1);

C = exp(B);
V_fit1 = (C./EPS_arr).*exp(A./EPS_arr);


plot(EPS_arr,V_fit1,'-g',LineWidth=2)
plot(EPS_arr,basin_volume,'.r',MarkerSize=20)
 %%  this is the same as before with out EPS^-1

% figure(30);
hold on 

delta = 1./EPS_arr(valume_index_cut:end); 

yy = log(basin_volume(valume_index_cut:end));

P2 = polyfit(delta,yy,1);

B = P2(2);
A = P2(1);

C = exp(B);
V_fit2 = (C.*exp(A./EPS_arr));

plot(EPS_arr,V_fit2,'-k',LineWidth=2)
% plot(EPS_arr,basin_volume,'.r',MarkerSize=20)
%% event function
function [check,stop,direction] = myeventfun(t,var,N)
check = var(N+1)>9; %prod(var(1:N) - 40);
stop = 1;  % Halt integration
direction = 0;
end
