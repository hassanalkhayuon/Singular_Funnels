% This code is created by Hassan Alkhayuon to the homoclinic bifurcation 
% in the eps vs alpha space
clear
warning off
addpath('..\')

% Colour-blind friendly colour
[red,yellow,green,blue ] = ...
    Colour_blind_friendly_colours();



alpha_init = 0.057;
EPSS_scan = linspace(0.057714914292180,8e-4,10);
for ind = 1:(length(EPSS_scan))
    EPSS = EPSS_scan(ind);
    gup_fun =@(alpha) homoclinic_gap(EPSS,alpha);
    alpha_init = fzero(gup_fun,alpha_init);
    alpha(ind) = alpha_init;
    ind
end

%%
alpha_1 = [alpha,1.0472];
EPSS_scan_1 = [EPSS_scan,0];

figure(2);
hold on
plot(alpha(1:end-1),EPSS_scan(1:end-1))
%%
% figure(1);
% cla
% homoclinic_plot(0.0005,0.96);

%% ODE function




