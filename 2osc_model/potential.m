% This code is created by Hassan Alkhayuon to
% calculate the potential for the avreaging system for mu in the 2 pahse
% oscilator system
% This is a Research project with Serhiy Yanchuk,
% Hildeberto Jardón-Kojakhmetov and Sebastian Wieczorek

function [V_mu] = potential(mu,par)
ome = par(1:2);
kappa = par(3);
eta = par(4);
alpha = par(5);

par = [ome;  kappa; eta; alpha; mu];

res_pot = 200;

mu_arr = linspace(-1,mu,res_pot);
g_bar_arr = zeros(size(mu_arr));

parfor ind = 1:res_pot
    g_bar_arr(ind) = mu_avg_2osc(mu_arr(ind),par);
end

V_mu = -1.*trapz(mu_arr,g_bar_arr); 

end