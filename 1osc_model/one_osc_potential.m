% This code is created by Hassan Alkhayuon to
% the potential of the avraged system for the osce phase oscelator 
% A Research project with Serhiy Yanchuk, Hildeberto Jardón-Kojakhmetov,
% and Sebastian Wieczorek

% the write hand side of the avraged system
eta = 10;
ome = -4;
alpha = pi/2;

% 
% figure(2);
% cla
% hold on 
% fplot(V1_,[3 5],'r')
% fplot(V2_,[-1 3],'-r')
% fplot(V2_,[5 15],'-r')
% plot([-1 15],[0 0],':k')

%% numarically 
f = @(m)f_avg(m,ome,alpha,eta);

res = 1000;

mu_arr = linspace(2, 16, res);

for ind = 2:res
    V(ind) = -1.*trapz( mu_arr(1:ind), f(mu_arr(1:ind)) );
end


figure(2);
hold on 
plot(...
    mu_arr,V,'-k')
%% avraging dynmaics 

function output = f_avg(m, ome, alpha, eta)

% Initialize the output array
output = zeros(size(m));

left_lim = -1 - ome;
right_lim = 1 - ome;


% Calculate output for each element of m
for ind_fun = 1:numel(m)
    if (m(ind_fun) < left_lim) || (m(ind_fun) > right_lim)
        output(ind_fun) =...
            -m(ind_fun) + eta .* (1 - (m(ind_fun) + ome) .* cos(alpha) - ...
            cos(alpha) .* sqrt((m(ind_fun) + ome).^2 - 1));
    else 
        output(ind_fun) = ...
            -m(ind_fun) + eta .* (1 - (m(ind_fun) + ome) .* cos(alpha) - ...
            sin(alpha) .* sqrt(1 - (m(ind_fun) + ome).^2));
    end
end
end