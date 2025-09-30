x_points = [3 4 7 10];
y_points = [-31 -29.5 -31 -42.6];

xx = linspace(x_points(1), x_points(end));
yy = interp1(x_points,y_points,xx,'spline');

% figure(20)
plot(xx,yy)