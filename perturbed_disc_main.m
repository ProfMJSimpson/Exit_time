%% Setup

clear; clc; close all
gmsh_path = 'C:\Users\simpsom3\gmsh-4.7.1-Windows64\gmsh.exe ';
Rc = 1; % radius of circle without perturbation
epsilon = 0.05; % perturbation parameter
g = @(theta) sin(3*theta) + cos(5*theta) - sin(theta); 
syms G(theta)
G(theta) = sin(3*theta);
R = @(theta) Rc*(1 + epsilon * g(theta)); % radius of peturbed geometry: [Equation (7)]
ref = 0.08; % refinement parameter in meshing (controls number of starting positions)
Nbpts = 1000; % number of points used to capture boundary in gmsh
P = 1; % probability of jumping at each step
Delta = 0.01; % step distance
tau = 1;  % step duration
D = P * Delta.^2 / (4*tau); % diffusivity
theta = linspace(0,2*pi,Nbpts);
N = 1000;
perturb_degree = 2;
fourier_terms = 25;

%% Mesh generation

perturbed_disc_mesh(theta,R,gmsh_path,ref,Nbpts);
mesh;
nodes = msh.POS(:,1:2);
elements = msh.TRIANGLES(:,1:3);
boundary_nodes = msh.LINES(find(msh.LINES(:,3)==99),1:2);
boundary_nodes = unique(reshape(boundary_nodes,2*size(boundary_nodes,1),1));

%% Random Walk Data

MeanExitTimes_Walk = perturbed_disc_walk(nodes(:, 1), nodes(:, 2), Delta, tau, P, N, R);

%% Finite Volume Data

MeanExitTimes_FVM = perturbed_disc_fvm(D, nodes, elements, boundary_nodes);

%% Perturbation Solution

% The function below also works for symbolic Rc, D, epsilon

T = perturbed_disc_perturbation(Rc, D, perturb_degree, fourier_terms, g, epsilon);

syms r theta
MeanExitTimes_Perturb = zeros(1, length(nodes(:, 1)));
for i =1:length(nodes(:, 1))
    if ~ismember(i, boundary_nodes)
        [ttheta, rr] = cart2pol(nodes(i,1),nodes(i,2));
        MeanExitTimes_Perturb(i) = double(subs(T, [r theta], [rr ttheta]));
    end
end

%% Tiled Layout: Figure 3
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
text(-1.67, -1.63, '(a)', 'FontSize', 28, 'Interpreter', 'latex')
t=subplot(1, 3, 2);
trisurf(elements, nodes(:, 1), nodes(:, 2), MeanExitTimes_FVM);
axis image
axis_specs(bottom, top);
set(t.Colorbar, 'Visible', 'off')
campos([0, 0, max(MeanExitTimes_FVM)]);
text(-1.67, -1.63, '(b)', 'FontSize', 28, 'Interpreter', 'latex')
t=subplot(1, 3, 3);
trisurf(elements, nodes(:, 1), nodes(:, 2), MeanExitTimes_Perturb)
axis image
axis_specs(bottom, top);
campos([0, 0, max(MeanExitTimes_Perturb)]);
text(-1.67, -1.63, '(c)', 'FontSize', 28, 'Interpreter', 'latex')

%% Detailed Perturbation Comparison 
bottom=0;top=10;
perturb_degree_1 = 1;   fourier_terms_1 = 1;
perturb_degree_2 = 1;   fourier_terms_2 = 10;
perturb_degree_3 = 2;   fourier_terms_3 = 10;
T1 = perturbed_disc_perturbation(Rc, D, perturb_degree_1, fourier_terms_1, g, epsilon);
T2 = perturbed_disc_perturbation(Rc, D, perturb_degree_2, fourier_terms_2, g, epsilon);
T3 = perturbed_disc_perturbation(Rc, D, perturb_degree_3, fourier_terms_3, g, epsilon);
%
syms r theta
MeanExitTimes_Perturb_1 = zeros(1, length(nodes(:, 1)));
for i =1:length(nodes(:, 1))
    if ~ismember(i, boundary_nodes)
        [ttheta, rr] = cart2pol(nodes(i,1),nodes(i,2));
        MeanExitTimes_Perturb_1(i) = double(subs(T1, [r theta], [rr ttheta]));
    end
end
MeanExitTimes_Perturb_2 = zeros(1, length(nodes(:, 1)));
for i =1:length(nodes(:, 1))
    if ~ismember(i, boundary_nodes)
        [ttheta, rr] = cart2pol(nodes(i,1),nodes(i,2));
        MeanExitTimes_Perturb_2(i) = double(subs(T2, [r theta], [rr ttheta]));
    end
end
MeanExitTimes_Perturb_3 = zeros(1, length(nodes(:, 1)));
for i =1:length(nodes(:, 1))
    if ~ismember(i, boundary_nodes)
        [ttheta, rr] = cart2pol(nodes(i,1),nodes(i,2));
        MeanExitTimes_Perturb_3(i) = double(subs(T3, [r theta], [rr ttheta]));
    end
end
Perturb_Error_1 = 100*abs(MeanExitTimes_FVM - MeanExitTimes_Perturb_1')/max(MeanExitTimes_FVM); % [Eq. 37]
Perturb_Error_2 = 100*abs(MeanExitTimes_FVM - MeanExitTimes_Perturb_2')/max(MeanExitTimes_FVM); % [Eq. 37]
Perturb_Error_3 = 100*abs(MeanExitTimes_FVM - MeanExitTimes_Perturb_3')/max(MeanExitTimes_FVM); % [Eq. 37]
%% Plotting: Figure 8
figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
t=subplot(1, 3, 1);
trisurf(elements, nodes(:, 1), nodes(:, 2), Perturb_Error_1);
axis image
axis_specs(bottom, top,true);
campos([0, 0, max(Perturb_Error_1)]);
set(t.Colorbar, 'Visible', 'off');
text(-1.67, -1.63, '(a)', 'FontSize', 28, 'Interpreter', 'latex')
t=subplot(1, 3, 2);
trisurf(elements, nodes(:, 1), nodes(:, 2), Perturb_Error_2);
axis image
axis_specs(bottom, top,true);
set(t.Colorbar, 'Visible', 'off')
campos([0, 0, max(Perturb_Error_2)]);
text(-1.67, -1.63, '(b)', 'FontSize', 28, 'Interpreter', 'latex')
t=subplot(1, 3, 3);
trisurf(elements, nodes(:, 1), nodes(:, 2), Perturb_Error_3)
axis image
axis_specs(bottom, top,true);
campos([0, 0, max(Perturb_Error_3)]);
text(-1.67, -1.63, '(c)', 'FontSize', 28, 'Interpreter', 'latex')

%%
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
    cb.YTick = linspace(bottom,top,5);
else
    title(cb, '$e(x, y)$', 'Interpreter', 'latex')
end
cb.Label.Interpreter='latex';
cb.TickLabelInterpreter='latex';
box on
grid off
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$y$', 'Interpreter', 'latex')
set(gcf, 'color', 'white') % for export_fig
set(gca, 'fontsize', 28, 'TickLabelInterpreter', 'latex')
end