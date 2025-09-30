
clear 
% System parameters

ome = [-3; -4];

eta = 10; % adaptive parameters
alpha = pi/2; % phase shift.
kappa = 1;


par = [ome;  kappa; eta; alpha];

plot_res = 200;

mu_arr = linspace(0, 15,plot_res);
V_mu = zeros(size(mu_arr));
tic
for ind = 1:plot_res
    % g_bar(ind) = mu_avg_2osc(mu_arr(ind),par);
    V_mu(ind) = potential(mu_arr(ind),par);
    disp(plot_res - ind)
end
toc
% figure(8);
% hold on
% plot(mu_arr,g_bar)
figure(9)
plot(mu_arr,V_mu)
xlim([mu_arr(1),mu_arr(end)])
ylim([-50 -20])
