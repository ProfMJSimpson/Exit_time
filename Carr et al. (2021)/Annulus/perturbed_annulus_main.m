%% Setup
clear; close all; clc

% GMSH Path
gmsh_path = '/Volumes/gmsh-4.7.1-MacOSX/Gmsh.app/Contents/MacOS/gmsh';
% gmsh_path = 'C:\Users\simpsom3\gmsh-4.7.1-Windows64\gmsh.exe';
if ~isfile(gmsh_path)
    warning('Download GMSH 4.7.1 from https://gmsh.info/bin/ and provide path to executable "gmsh".');
end

% Case = 'U1'; perturbed = false; perturb_degree = 0; % Figure 2
% Case = 'U2'; perturbed = false; perturb_degree = 0; % Figure 3
% Case = 'P1'; perturbed = true; perturb_degree = 2; % Figure 4
Case = 'P2'; perturbed = true; perturb_degree = 2; % Figure 5

% Symbolic variables
syms g1(theta) g2(theta) r

% Geometry parameters
R1 = 2; % unperturbed inner radius
R2 = 3; % unperturbed outer radius
epsilon = 0.05; % perturbation parameter
g1(theta) = sin(3*theta) + cos(5*theta); % perturbation function on inner radius
g2(theta) = cos(3*theta); % perturbation function on outer radius
if strcmp(Case,'U1')
    inner_type = 'absorb'; % inner boundary condition
    g1(theta) = 0; % perturbation function on inner radius
    g2(theta) = 0; % perturbation function on outer radius
elseif strcmp(Case,'U2')
    inner_type = 'reflect'; % inner boundary condition
    g1(theta) = 0; % perturbation function on inner radius
    g2(theta) = 0; % perturbation function on outer radius
elseif strcmp(Case,'P1')
    inner_type = 'absorb'; % inner boundary condition
    g1(theta) = sin(3*theta) + cos(5*theta); % perturbation function on inner radius
    g2(theta) = cos(3*theta); % perturbation function on outer radius
elseif strcmp(Case,'P2')
    inner_type = 'reflect'; % inner boundary condition
    g1(theta) = sin(3*theta) + cos(5*theta); % perturbation function on inner radius
    g2(theta) = cos(3*theta); % perturbation function on outer radius
end

% Random Walk Parameters
delta = 0.05; % step distance       
tau = 1;  % step duration
P = 1; % probability of jumping at each step
D = P * delta^2/(4*tau); % diffusivity
num_sims = 50; % number of random walk simulations

% FVM solution/meshing parameters
mesh_ref = 0.2; % characteristic length in gmsh
num_bpts = 1000; % number of points used to capture each boundary in gmsh
theta_vals = linspace(0, 2*pi, num_bpts); % generate initial angles

% Numerical parameters
fourier_terms = 25; % where to truncate Fourier series
% perturb_degree = 3; % degree of perturbation expansion

%% Mesh Generation
G_inner = matlabFunction(g1); % function handle rather than symbolic 
G_outer = matlabFunction(g2); % function handle rather than symbolic
R_inner = @(theta) R1 * (1 + epsilon * G_inner(theta)); % inner boundary
R_outer = @(theta) R2 * (1 + epsilon * G_outer(theta)); % outer boundary
perturbed_annulus_mesh(theta_vals,R_inner, R_outer,gmsh_path,mesh_ref,num_bpts);
Annulus;
nodes = msh.POS(:,1:2);
elements = msh.TRIANGLES(:,1:3);
outer_boundary_nodes = msh.LINES(find(msh.LINES(:,3)==55),1:2);
inner_boundary_nodes = msh.LINES(find(msh.LINES(:,3)==75),1:2);
outer_boundary_nodes = unique(reshape(outer_boundary_nodes,2*size(outer_boundary_nodes,1),1));
inner_boundary_nodes = unique(reshape(inner_boundary_nodes,2*size(inner_boundary_nodes,1),1));
boundary_nodes = union(outer_boundary_nodes, inner_boundary_nodes);

%% Plot Geometry
figure
triplot(elements, nodes(:, 1), nodes(:, 2))
box on
hold on
tint = linspace(0,2*pi,100); rint1 = R_inner(tint); rint2 = R_outer(tint);
[xint1,yint1] = pol2cart(tint,rint1);
plot(xint1,yint1,'r','LineWidth',2)
[xint2,yint2] = pol2cart(tint,rint2);
plot(xint2,yint2,'r','LineWidth',2)
axis image
drawnow

%% Random Walk
ExitTime_Walk = perturbed_annulus_walk(nodes(:, 1), nodes(:, 2), ...
    delta, tau, P, num_sims, R_inner, R_outer, inner_type);
disp('Random Walk complete')

%% FVM Solution
ExitTime_FVM = perturbed_annulus_fvm(D, nodes, elements, ...
    outer_boundary_nodes, inner_boundary_nodes, inner_type);
disp('FVM complete')

%% Perturbation Solution
if perturbed
    a = R1; b = R2;
    Perturbation = perturbed_annulus_perturbation(R1, R2, g1, g2, ...
        D, perturb_degree, fourier_terms, epsilon, inner_type, a, b, D);
    ExitTime_Perturbation = zeros(1, length(nodes(:, 1)));
    if strcmp(Case,'P1')
        absorbing_nodes = boundary_nodes;
    elseif strcmp(Case,'P2')
        absorbing_nodes = outer_boundary_nodes;
    end
    for i = 1:length(nodes(:, 1))
        if ~ismember(i, absorbing_nodes)
            [ttheta, rr] = cart2pol(nodes(i, 1),nodes(i, 2));
            ExitTime_Perturbation(i) = double(subs(Perturbation, [r theta], [rr ttheta]));
        end
    end
else
    [ttheta, rr] = cart2pol(nodes(:, 1),nodes(:, 2));
    if strcmpi(inner_type, 'reflect')
        ExitTime_Perturbation = (R2^2 - rr.^2 + 2*R1^2 * log(rr/R2))/(4*D);
    elseif strcmpi(inner_type, 'absorb')
        ExitTime_Perturbation = 1/(4*D*log(R1/R2)) * (R1^2*log(rr/R2) + ...
            R2^2*log(R1./rr) + rr.^2*log(R2/R1));
    end
end
disp('Perturbation complete')
fprintf('%g nodes and %g triangular elements\n',length(nodes(:,1)),length(elements(:,1)))

%% Plots
bottom = 0;
top = max([ExitTime_FVM]);
top = round(top, -2);

figure;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 1]);
t=subplot(1, 3, 1);
trisurf(elements, nodes(:, 1), nodes(:, 2), ExitTime_Walk);
axis image
axis_specs(bottom, top);
campos([0, 0, max(ExitTime_Walk)]);
axis([min([nodes(:,1);-R2]),max([nodes(:,1),;R2]),min([nodes(:,2);-R2]),max([nodes(:,2);R2])])
set(t.Colorbar, 'Visible', 'off');
text(-4.73, -4.46, '(a)', 'FontSize', 28, 'Interpreter', 'latex')
t=subplot(1, 3, 2);
trisurf(elements, nodes(:, 1), nodes(:, 2), ExitTime_Perturbation)
axis image
axis_specs(bottom, top);
set(t.Colorbar, 'Visible', 'off')
campos([0, 0, max(ExitTime_Perturbation)]);
axis([min([nodes(:,1);-R2]),max([nodes(:,1),;R2]),min([nodes(:,2);-R2]),max([nodes(:,2);R2])])
text(-4.73, -4.46, '(b)', 'FontSize', 28, 'Interpreter', 'latex')
t=subplot(1, 3, 3);
trisurf(elements, nodes(:, 1), nodes(:, 2), ExitTime_FVM);
axis image
axis_specs(bottom, top);
campos([0, 0, max(ExitTime_FVM)]);
axis([min([nodes(:,1);-R2]),max([nodes(:,1),;R2]),min([nodes(:,2);-R2]),max([nodes(:,2);R2])])
text(-4.73, -4.46, '(c)', 'FontSize', 28, 'Interpreter', 'latex')

% Sub-function (plot specifications)
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
set(gca,'XTick',-3:3:3,'YTick',-3:3:3)
if bar
    cb = colorbar;
    title(cb, '$T$', 'Interpreter', 'latex')
    cb.Label.Interpreter='latex';
    cb.TickLabelInterpreter='latex';
    yticks = linspace(bottom,top,5);
    cb.YTick = yticks;
    cb.Limits = [cb.YTick(1),cb.YTick(end)];
end
box on
grid off
xlabel('$x$', 'Interpreter', 'latex')
ylabel('$y$', 'Interpreter', 'latex')
set(gcf, 'color', 'white') % for export_fig
set(gca, 'fontsize', 28, 'TickLabelInterpreter', 'latex')
end