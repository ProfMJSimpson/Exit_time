%% Setup
clear; clc; close all
gmsh_path = 'C:\Users\simpsom3\gmsh-4.7.1-Windows64\gmsh.exe ';
Rc = 1; % radius of disc without perturbation
R = @(theta) 1; % constant; disc
ref = 0.08; % refinement parameter in meshing (controls number of starting positions)
Nbpts = 1000; % number of points used to capture boundary in gmsh
P = 1; % probability of jumping at each step
Delta = 1/100; % step distance
tau = 1;  % step duration
D = P * Delta.^2 / (4*tau); % diffusivity
theta = linspace(0,2*pi,Nbpts);
N = 1000; % number of realisations

%% Mesh generation

unperturbed_disc_mesh(theta,R,gmsh_path,ref,Nbpts);
mesh;
nodes = msh.POS(:,1:2);
elements = msh.TRIANGLES(:,1:3);
boundary_nodes = msh.LINES(find(msh.LINES(:,3)==99),1:2);
boundary_nodes = unique(reshape(boundary_nodes,2*size(boundary_nodes,1),1));

%% Random Walk Data

MeanExitTimes_Walk = unperturbed_disc_walk(nodes(:, 1), nodes(:, 2), Rc, Delta, tau, P, N);

%% Finite Volume Data

MeanExitTimes_FVM = unperturbed_disc_fvm(D, nodes, elements, boundary_nodes);

%% Exact Data

exact_soln = @(r) (Rc^2-r.^2)./(4*D); % [Eq. 4]
r_vals = sqrt(nodes(:, 1).^2 + nodes(:, 2).^2);
MeanExitTimes_Exact = exact_soln(r_vals);

%% Tiled Layout: Figure 1

bottom = 0;
top = 10000;

figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
t=subplot(1, 3, 1);
trisurf(elements, nodes(:, 1), nodes(:, 2), MeanExitTimes_Walk);
axis image
axis_specs(bottom, top);
campos([0, 0, max(MeanExitTimes_Walk)]);
set(t.Colorbar, 'Visible', 'off');
text(-1.63, -1.46, '(a)', 'FontSize', 28, 'Interpreter', 'latex')
t=subplot(1, 3, 2);
trisurf(elements, nodes(:, 1), nodes(:, 2), MeanExitTimes_FVM);
axis image
axis_specs(bottom, top);
set(t.Colorbar, 'Visible', 'off')
campos([0, 0, max(MeanExitTimes_FVM)]);
text(-1.63, -1.46, '(b)', 'FontSize', 28, 'Interpreter', 'latex')
t=subplot(1, 3, 3);
trisurf(elements, nodes(:, 1), nodes(:, 2), MeanExitTimes_Exact)
axis image
axis_specs(bottom, top);
campos([0, 0, max(MeanExitTimes_Exact)]);
text(-1.63, -1.46, '(c)', 'FontSize', 28, 'Interpreter', 'latex')

%%
function axis_specs(bottom,top,bar)
arguments
    bottom = 0
    top = 100
    bar = true
end
shading interp
lighting phong
camtarget([0 0 0])
caxis manual
caxis([bottom top]);
if bar
    cb = colorbar;
    title(cb, '$T$', 'Interpreter', 'latex')
    cb.Label.Interpreter='latex';
    cb.TickLabelInterpreter='latex';
    cb.YTick = linspace(bottom,top,5);
end
box on
grid off
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$y$', 'Interpreter', 'latex')
set(gcf, 'color', 'white') % for export_fig
set(gca, 'fontsize', 28, 'TickLabelInterpreter', 'latex')
end