g = @(theta) sin(3*theta) + cos(5*theta) - sin(theta);
h = @(theta) cos(3*theta) + sin(5*theta) - cos(theta);
a = 2;
b = 1;
figure
tiledlayout('flow')
nexttile
ellipp(3, 1, g, h, 0.005)
nexttile
ellipp(3, 1, g, h, 0.05)
nexttile
ellipp(3, 1, g, h, 0.5)
nexttile
ellipp(3, 1, g, h, 5)
nexttile
ellipp(3, 1, g, h, 7)
nexttile
ellipp(5, 2, g, h, 0.000071)
nexttile
ellipp(3, 2, g, h, 3.3321)
nexttile
ellipp(3, 1, g, h, 0.05)
nexttile
ellipp(5, 0.5, g, h, 0.1)
nexttile
ellipp(3, 2, g, h, 0.7)
%%
g = @(theta) sin(3*theta)+cos(5*theta)-sin(theta);
figure
tiledlayout('flow')
circc(1, g, 1/20)
nexttile
circc(1, g, 0.00005)
nexttile
circc(1, g, 3)
nexttile
circc(1, g, 20)
nexttile
circc(0.2, g, 0.01)
nexttile
circc(5, g, 5)
nexttile
circc(5, g, 6)
nexttile
circc(1, g, 6)
nexttile
circc(1, g, 0.3)
nexttile
circc(1, g, 0.729)


function ellipp(a, b, g, h, epsilon)
c1 = [sqrt(a^2-b^2) 0];
c2 = -c1;
sum_foci = 2*a;
theta = linspace(0,2*pi,1000);
x_pts = a*(1 + epsilon*g(theta)).*cos(theta);
y_pts = b*(1 + epsilon*h(theta)).*sin(theta);
points = [x_pts' y_pts'];
dist1 = (points(:, 1) - c1(1)).^2 + (points(:, 2) - c1(2)).^2;
dist2 = (points(:, 1) - c2(1)).^2 + (points(:, 2) - c2(2)).^2;
sum_foci_2 = sqrt(dist1) + sqrt(dist2);
max_diff = norm(sum_foci - sum_foci_2, Inf)
plot(points(:,1),points(:,2))
hold on
fimplicit(@(x, y) x^2/a^2+y^2/b^2-1)
title(sprintf('(a, b, $\\epsilon$, T) = (%g, %g, %g, %g)', a, b, epsilon, max_diff), 'Interpreter', 'latex')
axis image
end
function circc(R, g, epsilon)
Rfnc = @(theta) R*(1 + epsilon*g(theta)) - R;
fplot(Rfnc, [0 2*pi])
hold on
fplot(@(theta) 0, [0,2*pi])
vals = linspace(0,2*pi,1000);
vals_2 = Rfnc(vals);
max_err = norm(vals_2, Inf)
idx = find(abs(vals_2)==max_err);
plot([vals(idx) vals(idx)], [0 Rfnc(vals(idx))])
title(sprintf('(R, $\\epsilon$, T) = (%g, %g, %g)', R, epsilon, max_err), 'Interpreter', 'latex')
end