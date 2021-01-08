%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FITTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all; clc; clear
I = imread('taiwan_map.png');
% Image source (modified):
% https://d-maps.com/carte.php?num_car=5616&lang=en
% The modifications made are to provide some simplification to the boundary
% detail, such as removing the area below the Shuangliu National Forest
% Recreation Area and some rivers.

%% Boundary Fitting
boundary_spacing_parameter = 40;

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
x = x(1:spacing:end);
y = y(1:spacing:end);

%% Shift and Rotate Appropriately

% Compute Centroid
x_bar = mean(x);
y_bar = mean(y);

% Shift Boundary by Centroid
x = x - x_bar;
y = y - y_bar;

% Rotate Boundary so that a > b
boundary_data = [x y];
rotation_angle = pi/3; % counterclockwise, this angle is obtained by eye
rotation_matrix = [cos(rotation_angle) -sin(rotation_angle); ...
    sin(rotation_angle) cos(rotation_angle)];
boundary_data = rotation_matrix*boundary_data';
boundary_data = boundary_data';
x = boundary_data(:, 1);
y = boundary_data(:, 2);

%% Scale
new_pts = normalize_data_isotropically([x y]);
x = new_pts(:, 1);
y = new_pts(:, 2);

%% Construct Initial Ellipse
theta_vals = linspace(0, 2*pi, 1000);

% Develop Rotation and Translate
A_sampson = [x y];
[theta_fastguaranteed] = fast_guaranteed_ellipse_estimate(A_sampson);
geometricEllipseParameters = fromAlgebraicToGeometricParameters(theta_fastguaranteed);
% ^ [major axis, minor axis, x center, y center, orientation]
% Need to first translate center to the origin.
x_bar = geometricEllipseParameters(3);
y_bar = geometricEllipseParameters(4);
A_sampson(:, 1) = A_sampson(:, 1) - x_bar;
A_sampson(:, 2) = A_sampson(:, 2) - y_bar;
% Need to rotate back to be aligned with the coordinate axes.
phi = geometricEllipseParameters(end)-pi; 
rotated_A_sampson = [cos(phi) sin(phi); -sin(phi) cos(phi)] * A_sampson';
rotated_A_sampson = rotated_A_sampson';

% Recreate with New Transformation
sampson_x = rotated_A_sampson(:, 1);
sampson_y = rotated_A_sampson(:, 2);
A_sampson = [sampson_x sampson_y];
theta_fastguaranteed = fast_guaranteed_ellipse_estimate(A_sampson);
geometricEllipseParameters = fromAlgebraicToGeometricParameters(theta_fastguaranteed);

% Data
a_sampson = geometricEllipseParameters(1);
b_sampson = geometricEllipseParameters(2);
sampson_ellipse_x = a_sampson*cos(theta_vals);
sampson_ellipse_y = b_sampson*sin(theta_vals);

%% Construct Perturbed Boundaries
G_sampson = 9;
H_sampson = 9;
epsilon = 1/10;
unperturbed_theta = atan2(sampson_y/b_sampson, sampson_x/a_sampson);
g_matrix = epsilon*[ones(length(sampson_x),1) cos((1:G_sampson)'.*unperturbed_theta(1:length(sampson_x))')' sin((1:G_sampson)'.*unperturbed_theta(1:length(sampson_x))')'];
h_matrix = epsilon*[ones(length(sampson_x),1) cos((1:H_sampson)'.*unperturbed_theta(1:length(sampson_x))')' sin((1:H_sampson)'.*unperturbed_theta(1:length(sampson_x))')'];
g_rhs = sampson_x./(a_sampson*cos(unperturbed_theta)) - 1; % [Eq. 40]
h_rhs = sampson_y./(b_sampson*sin(unperturbed_theta)) - 1; % [Eq. 41]
g_solution = g_matrix \ g_rhs; % [Eq. 40]
h_solution = h_matrix \ h_rhs; % [Eq. 41]
syms g(theta) h(theta)
g_sampson(theta) = dot(g_solution, [1 cos((1:G_sampson)*theta) sin((1:G_sampson)*theta)]); % [Eq. 38]
h_sampson(theta) = dot(h_solution, [1 cos((1:H_sampson)*theta) sin((1:H_sampson)*theta)]); % [Eq. 39]
sampson_perturbed_x = a_sampson*(1+epsilon*g_sampson(theta_vals)).*cos(theta_vals); % [Eq. 16]
sampson_perturbed_y = b_sampson*(1+epsilon*h_sampson(theta_vals)).*sin(theta_vals); % [Eq. 17]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIRST PASSAGE TIME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Configuration
ref = 0.08;
bottom = 0;
top = 10000;
gmsh_path = 'C:\Users\simpsom3\gmsh-4.7.1-Windows64\gmsh.exe ';
Nbpts = 1000;
step_size = 0.01;
theta_vals = linspace(0, 2*pi, Nbpts);
time_step = 1;
move_prob = 1;
num_sims = 1000;
D = move_prob*step_size^2/(4*time_step);
perturb_degree_sampson = 2;
fourier_terms_sampson = 25;
G_sampson = matlabFunction(g_sampson); % [Eq. 38]
H_sampson = matlabFunction(h_sampson); % [Eq. 39]
x_func_sampson = @(psi) a_sampson*(1+epsilon*G_sampson(psi)).*cos(psi); % [Eq. 16]
y_func_sampson = @(psi) b_sampson*(1+epsilon*H_sampson(psi)).*sin(psi); % [Eq. 17]

%% Mesh Generation
perturbed_ellipse_mesh(sampson_perturbed_x, sampson_perturbed_y, gmsh_path, ref, Nbpts);
Ellipse;
sampson_nodes = msh.POS(:, 1:2);
sampson_elements = msh.TRIANGLES(:, 1:3);
sampson_boundary_nodes = msh.LINES(find(msh.LINES(:, 3) == 99), 1:2);
sampson_boundary_nodes = unique(reshape(sampson_boundary_nodes, 2*size(sampson_boundary_nodes, 1), 1));

%% FVM
FVM_sampson = perturbed_ellipse_fvm(D, sampson_nodes, sampson_elements, sampson_boundary_nodes);

%% Perturbation
perturbation_solution_sampson = perturbed_ellipse_perturbation(a_sampson, b_sampson, g_sampson, h_sampson, D, perturb_degree_sampson, fourier_terms_sampson, epsilon); % [Eq. 18]
perturbation_function_sampson = matlabFunction(perturbation_solution_sampson); % [Eq. 18]
perturbation_sampson = perturbation_function_sampson(sampson_nodes(:, 1), sampson_nodes(:, 2));

%% Error
bottom_error = 0;
top_error = 20;

error_sampson = 100*abs(FVM_sampson - perturbation_sampson)/max(FVM_sampson); % [Eq. 37]

%% Plots: Figure 7
main_figure = figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);

theta = pi-pi/3;
R = [1 0; 0 -1] * [cos(theta) -sin(theta); sin(theta) cos(theta)]; % rotate and flip to make map correct orientation
for i = 1:size(sampson_nodes, 1)
    sampson_nodes(i, 1:2) = (R * sampson_nodes(i, 1:2)');
end

sampson_fvm_plot = subplot(1, 3, 1);
trisurf(sampson_elements, sampson_nodes(:, 1), sampson_nodes(:, 2), FVM_sampson);
axis image
passage_axis(bottom, top);
campos([0 0 max(FVM_sampson)]);
set(sampson_fvm_plot.Colorbar, 'Visible', 'off');
text(-1.9, -2.35, '(a)', 'FontSize', 28, 'Interpreter', 'latex')

sampson_perturbation_plot = subplot(1, 3, 2);
trisurf(sampson_elements, sampson_nodes(:, 1), sampson_nodes(:, 2), perturbation_sampson);
axis image
passage_axis(bottom, top);
campos([0 0 max(perturbation_sampson)]);
text(-1.9, -2.35, '(b)', 'FontSize', 28, 'Interpreter', 'latex')

sampson_error_plot = subplot(1, 3, 3);
trisurf(sampson_elements, sampson_nodes(:, 1), sampson_nodes(:, 2), error_sampson);
axis image
passage_axis(bottom_error, top_error, true, true);
campos([0 0 max(error_sampson)]);
text(-1.8, -2.35, '(c)', 'FontSize', 28, 'Interpreter', 'latex')

% % Comparing to Unperturbed Ellipse
% Theory: Consider an ellipse x^2/a^2 + y^2/b^2, a > b. Then the foci are
% (c, 0) and (-c, 0), where c = (a^2 - b^2)^(1/2), and the sum of the
% distances from the foci to any point on the ellipse is constant,
% equalling 2a. Thus, we compare this value 2a to the maximum sum from the
% foci on the perturbed boundary.
% a = a_sampson;
% b = b_sampson;
% c_x = (a^2 - b^2)^(1/2);
% c_y = 0;
% c1 = [c_x c_y];
% c2 = -c1;
% sum_foci = 2*a;
% x_c = x_func_sampson;
% y_c = y_func_sampson;
% points = [x_c(theta_vals)' y_c(theta_vals)'];
% dist1 = (points(:, 1) - c1(1)).^2 + (points(:, 2) - c1(2)).^2;
% dist2 = (points(:, 1) - c2(1)).^2 + (points(:, 2) - c2(2)).^2;
% sum_foci_2 = sqrt(dist1) + sqrt(dist2);
% max_diff = norm(sum_foci - sum_foci_2, Inf)
% figure
% plot(points(:,1),points(:,2))
% hold on
% fimplicit(@(x, y) x^2/a^2+y^2/b^2-1)
% axis image
%% Axis Aesthetic Function
function passage_axis(bottom, top, err, yposlab)
arguments
    bottom = 0
    top = 100
    err = false
    yposlab = false
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
if ~yposlab
    ylabel('$y$', 'Interpreter', 'latex');
else
    text(-1.7, -0.3063, '$y$', 'Interpreter', 'latex', 'Rotation', 90, 'FontSize', 28)
end
set(gcf, 'color', 'white') % for export_fig
if err
    cb.YTick = 0:5:top;
end
set(gca, 'fontsize', 28, 'TickLabelInterpreter', 'latex')
end