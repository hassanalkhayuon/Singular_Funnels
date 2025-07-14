% This code is created by Hassan Alkhayuon to plot the Eps vs Alpha
% bifurcation diagram. This is the singular basin, a collaburation
% This is a Research project with Serhiy Yanchuk 
% and Hildeberto Jardón-Kojakhmetov

% load data 
load('../Data/BD_eps_vs_alpha_v2.mat')

% bifurcation curves
% eps_h1 = @(alp)interp1(alpha_arr_h1,EPS_arr_h1,alp,'spline');
eps_h1 = linspace(EPS_arr_h1(1),EPS_arr_h1(end),100);
eps_h2 = @(alp)interp1(alpha_arr_h2,EPS_arr_h2,alp,'spline');
eps_h3 = @(alp)interp1(alpha_arr_h3,EPS_arr_h3,alp,'spline');
eps_H1 = @(alp)interp1(alpha_arr_H1,EPS_arr_H1,alp,'spline');
eps_H2 = @(alp)interp1(alpha_arr_H2,EPS_arr_H2,alp,'spline');

% alpha_h1 = linspace(alpha_arr_h1(1),alpha_arr_h1(end),100);
alpha_h1 = @(EPSSS)interp1(EPS_arr_h1,alpha_arr_h1,EPSSS,'spline');
alpha_h2 = linspace(alpha_arr_h2(1),alpha_arr_h2(end),100);
alpha_h3 = linspace(alpha_arr_h3(1),alpha_arr_h3(end),100);
alpha_H1 = linspace(alpha_arr_H1(1),alpha_arr_H1(end),100);
alpha_H2 = linspace(alpha_arr_H2(1),alpha_arr_H2(end),100);

%% shaded regions
res = 1000;
alph_scan = linspace(-1.5,3,res);
eps_scan = linspace(0,0.1,res);

color_mat = ones(res,res);
for ind_alpha = 1:res
    alpha = alph_scan(ind_alpha);
    for ind_eps = 1:res
        EPS = eps_scan(ind_eps);
        % region 1
        if ((alpha<alpha_h2(end))&&(EPS>eps_h2(alpha)))
            color_mat(ind_eps,ind_alpha) = 0;
        end
        if (alpha>=alpha_h2(end))&&(alpha<alpha_h3(end))
            color_mat(ind_eps,ind_alpha) = 0;
        end
        if ( alpha>=alpha_h3(end) && alpha<alpha_h3(1) && ...
                EPS>eps_H2(alpha) )
            color_mat(ind_eps,ind_alpha) = 0;
        end
        if ( alpha>=alpha_h3(end) && EPS<eps_H2(alpha) && ...
                EPS>eps_h3(alpha) )
            color_mat(ind_eps,ind_alpha) = 0;
        end
        if ((alpha<alpha_h2(end))&&(EPS>eps_h2(alpha)))
            color_mat(ind_eps,ind_alpha) = 0;
        end
        if (alpha <= alpha_h1(EPS) )
            color_mat(ind_eps,ind_alpha) = 1;
        end
    end
end

%% plot
figure(3);
cla

pp = pcolor(alph_scan,eps_scan,color_mat);
pp.LineStyle = 'none';
pp.FaceAlpha = 0.1;
colormap([1 1 1;1 0 0]);
hold on 
plot(...
    alpha_h1(eps_h1), eps_h1,'-k','LineWidth',2)
plot(...
    alpha_h2,eps_h2(alpha_h2),'-k','LineWidth',2)
plot(...
    alpha_h3,eps_h3(alpha_h3),'-k','LineWidth',2)
plot(...
    alpha_H1,eps_H1(alpha_H1), '-b','LineWidth',2)
plot(...
    alpha_H2,eps_H2(alpha_H2), '-b','LineWidth',2)

plot(...
    alpha_h2(end),eps_h2(alpha_h2(end)),'.k',...
    alpha_h3(end),eps_h3(alpha_h3(end)),'.k',...
    'MarkerSize',30)

axis([-1.5 2.8 0 0.1])


