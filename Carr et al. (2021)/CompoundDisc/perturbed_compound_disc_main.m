%% Setup
clear; close all; clc

% GMSH Path
gmsh_path = '/Volumes/gmsh-4.7.1-MacOSX/Gmsh.app/Contents/MacOS/gmsh';
% gmsh_path = 'C:\Users\simpsom3\gmsh-4.7.1-Windows64\gmsh.exe';
if ~isfile(gmsh_path)
    warning('Download GMSH 4.7.1 from https://gmsh.info/bin/ and provide path to executable "gmsh".');
end

% Case = 'U1'; perturbed = false; perturb_degree = 0; % Figure 6
% Case = 'U2'; perturbed = false; perturb_degree = 0; % Figure 7
% Case = 'P1'; perturbed = true; perturb_degree = 2; % Figure 8
Case = 'P2'; perturbed = true; perturb_degree = 2; % Figure 9

% Geometry parameters
R1 = 2; % unperturbed inner radius 
R2 = 3; % unperturbed outer radius
epsilon = 0.05; % perturbation parameter
if strcmp(Case,'U1') || strcmp(Case,'U2')
    Gfunc = @(theta) 0; % perturbation function on interface
elseif strcmp(Case,'P1') || strcmp(Case,'P2')
    Gfunc = @(theta) sin(3*theta) + cos(5*theta); % perturbation function on interface
end

% Random walk parameters
delta = 0.05; % step distance
tau = 1; % step duration
if strcmp(Case,'U1') || strcmp(Case,'P1')
    P = [0.1, 1]; % probabilities of stepping in each layer
elseif strcmp(Case,'U2') || strcmp(Case,'P2')
    P = [1, 0.1]; % probabilities of stepping in each layer
end
D = P*delta^2/(4*tau); % diffusivities
n = 24; % number of random walk directions
num_sims = 50; % number of random walk simulations

% Perturbation solution parameters
fourier_terms = 25;

% FVM solution/meshing parameters
mesh_ref = 0.2;
num_bpts = 1000;

%% Meshing
R1star = @(theta) R1 * (1 + epsilon*Gfunc(theta));
R_vec = {@(theta) 0, R1star, @(theta) R2, @(theta) Inf};
gmsh = genmesh(R_vec, mesh_ref, num_bpts, gmsh_path);
nodes = gmsh.nodes;
elements = gmsh.elements;
bounds = gmsh.boundary_nodes;
x0 = nodes(:, 1);
y0 = nodes(:, 2);
[t0, r0] = cart2pol(x0, y0);
figure;
triplot(elements,x0,y0)
hold on
tint = linspace(0,2*pi,100); rint1 = R1star(tint); rint2 = R2*ones(size(tint));
[xint1,yint1] = pol2cart(tint,rint1);
plot(xint1,yint1,'r','LineWidth',2)
[xint2,yint2] = pol2cart(tint,rint2);
plot(xint2,yint2,'r','LineWidth',2)
axis image
drawnow

%% Random Walk
walks = zeros(1, size(gmsh.nodes, 1));
for i = 1:length(x0)
    if(~ismember(i, bounds))
        x0elem = x0(i);
        y0elem = y0(i);     
        for j = 1:num_sims
            %walks(i) = walks(i) + het_disc(R_vec, P, [x0elem, y0elem], delta, tau, n, false);
            walks(i) = walks(i) + het_disc2(R1star, R2, P, [x0elem, y0elem], delta, tau, n, false);
        end   
        walks(i) = walks(i) / num_sims;
    end
    clc
    disp(['Random Walk completion: ',num2str(round(100*i/length(x0))),'%']);
end
%% FVM Solution
fvmmesh = fvm_properties(gmsh, R_vec);
fvmsol = full(fvm_heterogeneous(D, fvmmesh));
disp('FVM complete')

%% Perturbation Solution
[tvals, rvals] = cart2pol(nodes(:,1), nodes(:,2));
if perturbed
    STfunc = pertcalc(R1, R2, D, Gfunc, 0, fourier_terms);
    if perturb_degree > 1
        STfunc = pertcalc(R1, R2, D, Gfunc, perturb_degree, fourier_terms, STfunc);
    end
    pertsols = evalperturb(x0, y0, R1star, STfunc, epsilon);
else
    pertsols = zeros(length(nodes(:,1)),1);
    for i = 1:length(nodes(:,1))
        if rvals(i) < R1
            pertsols(i) = (R1^2 - rvals(i).^2)/(4*D(1)) + (R2^2 - R1^2)/(4*D(2));
        else
            pertsols(i) = (R2^2 - rvals(i).^2)/(4*D(2));
        end
    end
end
disp('Perturbation complete')
fprintf('%g nodes and %g triangular elements\n',length(nodes(:,1)),length(elements(:,1)))

%% Plots
ExitTime_Walk = walks;
ExitTime_FVM = fvmsol;
ExitTime_Exact = pertsols(:,perturb_degree+1);
bottom = 0;
top = max([ExitTime_FVM; ExitTime_Exact]);
top = round(top, -2);
figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 1, 1]);
t=subplot(1, 3, 1);
trisurf(elements, nodes(:, 1), nodes(:, 2), ExitTime_Walk);
hold on
plot3(xint1,yint1,0.99999*max(ExitTime_Walk)*ones(size(xint1)),'k','LineWidth',2)
axis image
axis_specs(bottom, top);
axis([-R2 R2 -R2 R2])
set(t.Colorbar, 'Visible', 'off');
campos([0, 0, max(ExitTime_Walk)]);
text(-4.73, -4.46, '(a)', 'FontSize', 28, 'Interpreter', 'latex')
t=subplot(1, 3, 2);
trisurf(elements, nodes(:, 1), nodes(:, 2), ExitTime_Exact)
axis image
axis_specs(bottom, top);
axis([-R2 R2 -R2 R2])
set(t.Colorbar, 'Visible', 'off')
campos([0, 0, max(ExitTime_Exact)]);
text(-4.73, -4.46, '(b)', 'FontSize', 28, 'Interpreter', 'latex')
t=subplot(1, 3, 3);
trisurf(elements, nodes(:, 1), nodes(:, 2), ExitTime_FVM);
axis image
axis_specs(bottom, top);
axis([-R2 R2 -R2 R2])
campos([0, 0, max(ExitTime_FVM)]);
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
    title(cb, '$T^{(1)}$, $T^{(2)}$', 'Interpreter', 'latex')
    cb.Label.Interpreter='latex';
    cb.TickLabelInterpreter='latex';
    yticks = linspace(bottom,top,5);
    cb.YTickLabel = cellstr(num2str(reshape(yticks,[],1),'%g'));
%     for k = 1:5
%         yticks(k) = round(yticks(k),max(floor(log10(yticks(k)))-1,2),'significant');
%     end
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