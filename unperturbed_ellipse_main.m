%% Setup

clc; clear; close all
a = 2; % semi-major axis
b = 1; % semi-minor axis
gmsh_path = 'C:\Users\simpsom3\gmsh-4.7.1-Windows64\gmsh.exe ';
Nbpts = 1000;
ref = 0.08;
theta_vals = linspace(0, 2*pi, Nbpts);
x_bnd = a*cos(theta_vals);
y_bnd = b*sin(theta_vals);
step_size = 1/100;                                                     
time_step = 1;                                                   
move_prob = 1;                                                  
num_sims = 1000;   
D = move_prob * step_size^2 / (4*time_step);  

%% Mesh Generation

unperturbed_ellipse_mesh(x_bnd, y_bnd, gmsh_path, ref, Nbpts)
Ellipse;
nodes = msh.POS(:,1:2);
elements = msh.TRIANGLES(:,1:3);
boundary_nodes = msh.LINES(find(msh.LINES(:,3)==99),1:2);
boundary_nodes = unique(reshape(boundary_nodes,2*size(boundary_nodes,1),1));

%% Exact Solution 

T = @(x, y) a^2*b^2/(2*D*(a^2+b^2))*(1-x.^2/a^2-y.^2/b^2); % [Equation (6)]
MeanExitTimes_Exact = T(nodes(:, 1), nodes(:, 2));

%% FVM

MeanExitTimes_FVM = unperturbed_ellipse_fvm(D, nodes, elements, boundary_nodes);
 
%% Random Walk

MeanExitTimes_Walk = unperturbed_ellipse_walk(nodes(:, 1), nodes(:, 2), step_size, time_step, move_prob, num_sims, a, b);

%% Tiled Layout: Figure 2
bottom=0;
top=10000;
figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
t=subplot(1, 3, 1);
trisurf(elements, nodes(:, 1), nodes(:, 2), MeanExitTimes_Walk);
axis image
axis_specs(bottom, top);
campos([0, 0, max(MeanExitTimes_Walk)]);
set(t.Colorbar, 'Visible', 'off');
text(-3.35, -2, '(a)', 'FontSize', 28, 'Interpreter', 'latex')
t=subplot(1, 3, 2);
trisurf(elements, nodes(:, 1), nodes(:, 2), MeanExitTimes_FVM);
axis image
axis_specs(bottom, top);
set(t.Colorbar, 'Visible', 'off')
campos([0, 0, max(MeanExitTimes_FVM)]);
text(-3.35, -2, '(b)', 'FontSize', 28, 'Interpreter', 'latex')
t=subplot(1, 3, 3);
trisurf(elements, nodes(:, 1), nodes(:, 2), MeanExitTimes_Exact)
axis image
axis_specs(bottom, top);
campos([0, 0, max(MeanExitTimes_Exact)]);
text(-3.35, -2, '(c)', 'FontSize', 28, 'Interpreter', 'latex')

%% Axis Aesthetics
function axis_specs(bottom,top)
shading interp
lighting phong
camtarget([0 0 0])
caxis manual
caxis([bottom top]);
cb = colorbar;
title(cb, '$T$', 'Interpreter', 'latex')
cb.Label.Interpreter='latex';
cb.TickLabelInterpreter='latex';
cb.YTick = linspace(bottom,top,5);
box on
grid off
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$y$', 'Interpreter', 'latex')
set(gcf, 'color', 'white') % for export_fig
set(gca, 'fontsize', 28, 'TickLabelInterpreter', 'latex')
end