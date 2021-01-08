%% Parameters

clear; clc; close all
syms g(theta) h(theta) n x_fnc(psi) y_fnc(psi)
a = 2;
b = 1;
gmsh_path = 'C:\Users\simpsom3\gmsh-4.7.1-Windows64\gmsh.exe ';
Nbpts = 1000;
ref = 0.08; 
epsilon = 0.05;
g(theta) = sin(3*theta) + cos(5*theta) - sin(theta);
h(theta) = cos(3*theta) + sin(5*theta) - cos(theta);
theta_vals = linspace(0, 2*pi, Nbpts);
x_bnd = a*(1 + epsilon * g(theta_vals)).* cos(theta_vals); % [Eq. 16]
y_bnd = b*(1 + epsilon * h(theta_vals)).* sin(theta_vals); % [Eq. 17]
G = matlabFunction(g); % [Eq. 16]
H = matlabFunction(h); % [Eq. 17]
x_fnc = @(psi) a*(1+epsilon*G(psi)).*cos(psi); % [Eq. 16]
y_fnc = @(psi) b*(1+epsilon*H(psi)).*sin(psi); % [Eq. 17]
x_bnd_1 = a*cos(theta_vals); 
y_bnd_1 = b*sin(theta_vals);
step_size = 0.01;                                                     
time_step = 1;                                                   
move_prob = 1;                                                  
num_sims = 1000;   
perturb_degree = 2;
N = 25;
D = move_prob * step_size^2 / (4*time_step);  

%% Mesh Generation

perturbed_ellipse_mesh(x_bnd, y_bnd, gmsh_path, ref, Nbpts)
Ellipse;
nodes = msh.POS(:,1:2);
elements = msh.TRIANGLES(:,1:3);
boundary_nodes = msh.LINES(find(msh.LINES(:,3)==99),1:2);
boundary_nodes = unique(reshape(boundary_nodes,2*size(boundary_nodes,1),1));


%% FVM

MeanExitTimes_FVM = perturbed_ellipse_fvm(D, nodes, elements, boundary_nodes);

%% Random Walk

MeanExitTimes_Walk = perturbed_ellipse_walk_2(nodes(:,1), nodes(:,2), step_size,time_step, move_prob, num_sims, nodes, boundary_nodes);

%% Perturbation

Perturbation_Soln = perturbed_ellipse_perturbation(a, b, g, h, D, perturb_degree, N, epsilon); % [Eq. 18]
T_perturb_func = matlabFunction(Perturbation_Soln); % [Eq. 18]
MeanExitTimes_Perturb = T_perturb_func(nodes(:, 1), nodes(:, 2));

%% Comparison: Comparing to FVM

Perturb_Error = 100*abs(MeanExitTimes_FVM - MeanExitTimes_Perturb)/max(MeanExitTimes_FVM); % [Eq. 37]

%% Tiled Layout (Figure 4)
bottom=0;
top= 10000;
figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
t=subplot(1, 3, 1);
trisurf(elements, nodes(:, 1), nodes(:, 2), MeanExitTimes_Walk);
axis image
axis_specs(bottom, top);
campos([0, 0, max(MeanExitTimes_Walk)]);
set(t.Colorbar, 'Visible', 'off');
text(-3.1, -2, '(a)', 'FontSize', 28, 'Interpreter', 'latex')
t=subplot(1, 3, 2);
trisurf(elements, nodes(:, 1), nodes(:, 2), MeanExitTimes_FVM);
axis image
axis_specs(bottom, top);
set(t.Colorbar, 'Visible', 'off')
campos([0, 0, max(MeanExitTimes_FVM)]);
text(-3.1, -2, '(b)', 'FontSize', 28, 'Interpreter', 'latex')
t=subplot(1, 3, 3);
trisurf(elements, nodes(:, 1), nodes(:, 2), MeanExitTimes_Perturb)
axis image
axis_specs(bottom, top);
campos([0, 0, max(MeanExitTimes_Perturb)]);
text(-3.1, -2, '(c)', 'FontSize', 28, 'Interpreter', 'latex')
   
%% Detailed Perturbation Comparison
bottom=0;top=10;
perturb_degree_1 = 1;   fourier_terms_1 = 1;
perturb_degree_2 = 1;   fourier_terms_2 = 10;
perturb_degree_3 = 2;   fourier_terms_3 = 10;
T1 = perturbed_ellipse_perturbation(a, b, g, h, D, perturb_degree_1, fourier_terms_1, epsilon); % [Eq. 18]
T2 = perturbed_ellipse_perturbation(a, b, g, h, D, perturb_degree_2, fourier_terms_2, epsilon); % [Eq. 18]
T3 = perturbed_ellipse_perturbation(a, b, g, h, D, perturb_degree_3, fourier_terms_3, epsilon); % [Eq. 18]
T_perturb_func_1 = matlabFunction(T1); % [Eq. 18]
T_perturb_func_2 = matlabFunction(T2); % [Eq. 18]
T_perturb_func_3 = matlabFunction(T3); % [Eq. 18]
MeanExitTimes_Perturb_1 = T_perturb_func_1(nodes(:, 1), nodes(:, 2));
MeanExitTimes_Perturb_2 = T_perturb_func_2(nodes(:, 1), nodes(:, 2));
MeanExitTimes_Perturb_3 = T_perturb_func_3(nodes(:, 1), nodes(:, 2));
Perturb_Error_1 = 100*abs(MeanExitTimes_FVM - MeanExitTimes_Perturb_1)/max(MeanExitTimes_FVM); % [Eq. 37]
Perturb_Error_2 = 100*abs(MeanExitTimes_FVM - MeanExitTimes_Perturb_2)/max(MeanExitTimes_FVM); % [Eq. 37]
Perturb_Error_3 = 100*abs(MeanExitTimes_FVM - MeanExitTimes_Perturb_3)/max(MeanExitTimes_FVM); % [Eq. 37]
%% Plotting: Figure 9
figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
t=subplot(1, 3, 1);
trisurf(elements, nodes(:, 1), nodes(:, 2), Perturb_Error_1);
axis image
axis_specs(bottom, top,true);
campos([0, 0, max(Perturb_Error_1)]);
set(t.Colorbar, 'Visible', 'off');
text(-3.1, -2, '(a)', 'FontSize', 28, 'Interpreter', 'latex')
t=subplot(1, 3, 2);
trisurf(elements, nodes(:, 1), nodes(:, 2), Perturb_Error_2);
axis image
axis_specs(bottom, top,true);
set(t.Colorbar, 'Visible', 'off')
campos([0, 0, max(Perturb_Error_2)]);
text(-3.1, -2, '(b)', 'FontSize', 28, 'Interpreter', 'latex')
t=subplot(1, 3, 3);
trisurf(elements, nodes(:, 1), nodes(:, 2), Perturb_Error_3)
axis image
axis_specs(bottom, top,true);
campos([0, 0, max(Perturb_Error_3)]);
text(-3.1, -2, '(c)', 'FontSize', 28, 'Interpreter', 'latex')

%% Axis Aesthetics
function axis_specs(bottom,top,err)
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
    cb.YTick = 0:5:10;
end
set(gca, 'fontsize', 28, 'TickLabelInterpreter', 'latex')
end