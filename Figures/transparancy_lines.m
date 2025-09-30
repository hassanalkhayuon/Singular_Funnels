h = findobj(gca,'Type','Line');
for ind = 1:12
h(ind).Color = [h(ind).Color(1:3), 0.3;]
end