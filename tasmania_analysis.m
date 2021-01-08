%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FITTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clc; clear
I = imread('tasmania_map.png');
% Image source (modified):
% http://mapsof.net/tasmania/blank-map-of-tasmania
% The modifications made are to provide some simplification to the boundary
% detail. Some of these are:
%       - Removing regions near and including Stanley.
%       - Removing the part on the north-west section of Woolnorth.
%       - Slightly smoothing out the Macquarie Heads, connecting with
%       Strahan to fill in the Macquarie Harbour Historic Site.
%       - Filling in the region by the Southwest National Park and Huon
%       Valley.
%       - Removed the section off Dunalley.
%       - Removed the section off Lauderdale.
%       - Filled in the River Tamar.

%% Extract Boundary
boundary_spacing_parameter = 30; % larger => more points

% Split Channels
red_channel = I(:, :, 1);
green_channel = I(:, :, 2);
blue_channel = I(:, :, 3);

% Detect Edges in Each Channel
edge_detection_method = 'Sobel';
edges_red_channel = edge(red_channel, edge_detection_method);
edges_green_channel = edge(green_channel, edge_detection_method);
edges_blue_channel = edge(blue_channel, edge_detection_method);

% Detect Edges
edges = edges_red_channel | edges_green_channel | edges_blue_channel;
se = strel('disk', 3);
edges = imclose(edges, se);

% Detect Boundaries
[boundaries, boundary_labels] = bwboundaries(edges);

% Boundary
boundary = boundaries{2};
x = boundary(:, 2);
y = boundary(:, 1);

% Sample Less Points
spacing = ceil(length(x) / boundary_spacing_parameter);
xx = x;
yy = y;
x = x(1:spacing:end);
y = y(1:spacing:end);

%% Shift Appropriately

% Compute Centroid
x_bar = mean(x);
xx_bar = mean(xx);
y_bar = mean(y);
yy_bar = mean(yy);

% Shift Boundary by Centroid
x = x - x_bar;
xx = xx - xx_bar;
y = y - y_bar;
yy = yy - y_bar;

%% Scale Image Isotropically
new_pts = normalize_data_isotropically([x y]);
x = new_pts(:, 1);
y = new_pts(:, 2);

new_pts_2 = normalize_data_isotropically([xx yy]);
xx = new_pts_2(:, 1);
yy = new_pts_2(:, 2);

%% Construct Initial Disc I
theta_vals = linspace(0, 2*pi, 1000);
spline_values = cscvn([x'; y']);
spline_form = fn2fm(spline_values, 'B-');
spline_vals = fnval(spline_form, spline_form.knots);
spline_x = spline_vals(1, :)';
spline_y = spline_vals(2, :)';
system_size_spline = length(spline_x);
A_spline = [spline_x spline_y];
coefficient_matrix_spline = [A_spline ones(system_size_spline, 1)];
rhs_vector_spline = A_spline(:, 1).^2 + A_spline(:, 2).^2; 
spline_soln = coefficient_matrix_spline \ rhs_vector_spline;
A = spline_soln(1);
B = spline_soln(2);
spline_x_c = A/2;
spline_y_c = B/2;

%% Shift by Circle Center
spline_x = spline_x - spline_x_c;
spline_y = spline_y - spline_y_c;
new_spline_boundary_x = x - spline_x_c;
new_spline_boundary_y = y - spline_y_c;

%% Construct Initial Disc II
spline_values = cscvn([x'; y']);
spline_form = fn2fm(spline_values, 'B-');
spline_vals = fnval(spline_form, spline_form.knots);
system_size_spline = length(spline_x);
A_spline = [spline_x spline_y];
coefficient_matrix_spline = [A_spline ones(system_size_spline, 1)]; 
rhs_vector_spline = A_spline(:, 1).^2 + A_spline(:, 2).^2;
spline_soln = coefficient_matrix_spline \ rhs_vector_spline; 
A = spline_soln(1);
B = spline_soln(2);
C = spline_soln(3);
spline_x_c = A/2;
spline_y_c = B/2;
R_spline = sqrt(4*C+A^2+B^2)/2;
spline_disc_x = R_spline*cos(theta_vals); 
spline_disc_y = R_spline*sin(theta_vals);

%% Construct Perturbed_Boundaries
epsilon = 1/10;
G_spline = 3;
unperturbed_theta = atan2(spline_y, spline_x);
g_matrix = epsilon*[ones(system_size_spline, 1) cos((1:G_spline)'.*unperturbed_theta(1:system_size_spline)')' sin((1:G_spline)'.*unperturbed_theta(1:system_size_spline)')']; % [Eq. 35]
g_rhs = 1/R_spline*sqrt(spline_x.^2+spline_y.^2)-1; % [Eq. 35]
g_solution = g_matrix \ g_rhs; % [Eq. 36]
syms g_spline(theta)
g_spline(theta) = dot(g_solution, [1 cos((1:G_spline)*theta) sin((1:G_spline)*theta)]); % [Eq. 34]
spline_perturbed_x = R_spline*(1+epsilon*g_spline(theta_vals)).*cos(theta_vals); % [Eq. 7]
spline_perturbed_y = R_spline*(1+epsilon*g_spline(theta_vals)).*sin(theta_vals); % [Eq. 7]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIRST PASSAGE TIME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Configuration
ref = 0.08;
bottom = 0;
top = 15000;
gmsh_path = 'C:\Users\simpsom3\gmsh-4.7.1-Windows64\gmsh.exe ';
Nbpts = 1000;
P = 1;
Delta = 0.01;
tau = 1;
D = P*Delta^2/(4*tau);
theta_vals = linspace(0, 2*pi,Nbpts);
perturb_degree_spline = 2; 
fourier_terms_spline = 25;
G_spline = matlabFunction(g_spline);
R_func_spline = @(theta) R_spline*(1+epsilon*g_spline(theta)); % [Eq. 7]

%% Mesh Generation
perturbed_disc_mesh(theta_vals,R_func_spline,gmsh_path,ref,Nbpts);
mesh;
spline_nodes = msh.POS(:,1:2);
spline_elements = msh.TRIANGLES(:,1:3);
spline_boundary_nodes = msh.LINES(find(msh.LINES(:,3)==99),1:2);
spline_boundary_nodes = unique(reshape(spline_boundary_nodes,2*size(spline_boundary_nodes,1),1));

%% FVM
FVM_spline = perturbed_disc_fvm(D, spline_nodes, spline_elements, spline_boundary_nodes);

%% Perturbation
T_spline = perturbed_disc_perturbation(R_spline, D, perturb_degree_spline, fourier_terms_spline, g_spline, epsilon);
syms r theta
perturbation_spline = zeros(1, length(spline_nodes(:, 1)));
for i=1:length(spline_nodes(:, 1))
    if ~ismember(i, spline_boundary_nodes)
        [ttheta, rr] = cart2pol(spline_nodes(i, 1), spline_nodes(i, 2));
        perturbation_spline(i) = double(subs(T_spline, [r theta], [rr ttheta]));
    end
end

%% Error
bottom_error = 0;
top_error = 20;

error_spline = 100*abs(FVM_spline - perturbation_spline')/max(FVM_spline); % [Eq. 37]

%% Plots: Figure 6
main_figure = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

theta = -pi;
R = [-1 0; 0 1] * [cos(theta) -sin(theta); sin(theta) cos(theta)]; % rotate and flip to make map correct orientation
for i = 1:size(spline_nodes, 1)
    spline_nodes(i, 1:2) = (R * spline_nodes(i, 1:2)');
end

% FVM

spline_fvm_plot = subplot(1, 3, 1);
trisurf(spline_elements, spline_nodes(:, 1), spline_nodes(:, 2), FVM_spline);
axis image
passage_axis(bottom, top);
campos([0 0 max(FVM_spline)]);
set(spline_fvm_plot.Colorbar, 'Visible', 'off');
text(-2.3, -2.55, '(a)', 'FontSize', 28, 'Interpreter', 'latex')

% Perturbation

spline_perturbation_plot = subplot(1, 3, 2);
trisurf(spline_elements, spline_nodes(:, 1), spline_nodes(:, 2), perturbation_spline);
axis image
passage_axis(bottom, top);
campos([0 0 max(perturbation_spline)]);
text(-2.3, -2.55, '(b)', 'FontSize', 28, 'Interpreter', 'latex')

% Error

spline_error_plot = subplot(1, 3, 3);
trisurf(spline_elements, spline_nodes(:, 1), spline_nodes(:, 2), error_spline);
axis image
passage_axis(bottom_error, top_error, true);
campos([0 0 max(error_spline)]);
text(-2.3, -2.55, '(c)', 'FontSize', 28, 'Interpreter', 'latex')
% % Comparing to Unperturbed Disc
% R_func_spline = @(theta) R_spline*(1+G_spline(theta));
% r_vals = R_func_spline(linspace(0,2*pi,1000));
% max_diff = norm(r_vals - R_spline, Inf)
% [x,y] = pol2cart(linspace(0,2*pi,1000),r_vals);
% figure
% plot(x, y)
% hold on
% fimplicit(@(x, y) x^2+y^2-R_spline)
% axis image
%%
% figure
fn = @(theta) R_func_spline(theta) - R_spline;
% fplot(fn, [0 2*pi])
% hold on
% fplot(@(theta) 0,[0,2*pi])
vals = linspace(0,2*pi,1000);
vals_2 = fn(vals);
max_err = norm(vals_2, Inf));
vpa(max_err);
% idx=find(abs(vals_2) == max_err);
% plot([vals(idx) vals(idx)], [0 fn(vals(idx))])
%% Axis Aesthetic Function
function passage_axis(bottom, top, err)
arguments
    bottom = 0
    top = 100
    err = false
end
shading interp
lighting phong
camtarget([0 0 0])
caxis manual
caxis([bottom top]);
cb = colorbar;
if ~err
    title(cb, '$T$', 'Interpreter', 'latex')
else
    title(cb, '$e(x, y)$', 'Interpreter', 'latex')
end
cb.Label.Interpreter='latex';
cb.TickLabelInterpreter='latex';
cb.YTick = linspace(bottom, top, 5);
box on
grid off
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$y$', 'Interpreter', 'latex')
set(gcf, 'color', 'white') % for export_fig
if err
    cb.YTick = 0:5:top;
end
set(gca, 'fontsize', 28, 'TickLabelInterpreter', 'latex')
end