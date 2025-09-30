% This code is created by Hassan Alkhayuon to
% compute the size of basin of attraction for 2 coupled adaptive phase
% oscillators by monte carlo on varying Delta_mu and epsilon
% A Research project with Serhiy Yanchuk, Hildeberto Jardón-Kojakhmetov,
% and sebastian Wieczorek

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

Eps = 0.05;
tend = 150;

Del_arr = [0:0.2:0.8,0.9,1,1.1,1.2:0.2:2.4];
Eps_arr = linspace(1e-3,1e-1,100);

res_Del = length(Del_arr);
res_EPS = length(Eps_arr);

basin_volume_l = zeros(size(Del_arr));

number_of_simulations = 100000;

rr = rand(number_of_simulations, 3);

tic
for ind_Del = 1:res_Del
    ome1 = ome2 + Del_arr(ind_Del);
    for ind_Eps = 1:res_EPS

        basin_vol_l = 0;

        Eps = Eps_arr(ind_Eps);
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
            end
        end
        basin_volume_l(ind_Del,ind_Eps) = basin_vol_l/number_of_simulations;

        disp([res_Del - ind_Del, res_EPS - ind_Eps])
    end
    save("Basin_vul_2osc_scan.mat")
end
toc
%% plotting
figure(8);
clf
pp=pcolor(Eps_arr, Del_arr, basin_volume_l);
pp.LineStyle = 'none';

%% plotting one curve 
figure(9);

gif_filename = 'basin_vol_vary_del';
 
for ind_Del = 1:res_Del-1
    clf
    hold on
    plot(Eps_arr,basin_volume_l(1:end-1,:)','-','Color',[0.6 0.6 0.6])
    plot(Eps_arr,basin_volume_l(ind_Del,:)','-k','LineWidth',2)
    yscale log
    axis([1e-2 1e-1 1e-4 1e0])
    % legend(['$\Delta_\omega =$' num2str(Del_arr(ind_Del))],)
    str_label = ['$\Delta_\omega = $' num2str(Del_arr(ind_Del))];
    annotation('textbox',...
    [0.66 0.19 0.085 0.053],...
    'String', str_label,...
    'LineStyle','none',...
    'FontSize',15);
    set(gca,'FontSize',15);
    xlabel('$\varepsilon$')
    ylabel('Basin volume')
    box on
    drawnow()
    
    % --- Capture frame as image ---
    frame = getframe(gcf);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256);

    % --- Write to GIF ---
    if ind_Del == 1
        imwrite(A,map,gif_filename,"gif","LoopCount",Inf,"DelayTime",1);
    else
        imwrite(A,map,gif_filename,"gif","WriteMode","append","DelayTime",1);
    end
    
end
%%

colors = [
    1, 0, 0, 0.2;        % Red
    0, 0.5, 0, 0.2;      % Dark Green
    0, 0, 1, 0.2;        % Blue
    0, 0, 0, 1;        % Black     
    1, 0, 1, 0.2;        % Magenta
    0.8, 0.8, 0, 0.2;    % Darker Yellow
    0.18,0.75,0.94, 1;   % Cyan
    0.5, 0.7, 1, 0.2;    % Light Blue
    0.5, 0.5, 0.5, 0.2;  % Gray
    1, 0.5, 0, 0.2;      % Orange
    0.5, 0, 0.5, 0.2;    % Purple
    0.6, 0.3, 0, 0.2;    % Brown
    1, 0.75, 0.8, 0.2;   % Pink
    0.5, 0.5, 0, 1     % Olive
];

figure(11);
cla
hold on
for ind_Del = 9:14
    plot(Eps_arr,basin_volume_l(ind_Del,:)',color = colors(ind_Del,:),...
        LineWidth =0.5)
    % legendEntries{ind_Del} = ['$\Delta_\omega =$' num2str(Del_arr(ind_Del))];
end
yscale log
box on
% legend(legendEntries, 'Location', 'southeast', 'NumColumns', 2);

%% event function

function [check, stop, direction] = myeventfun(t, var, N)
check = var(N+1)>8;
stop = 1;
direction = 0;
end