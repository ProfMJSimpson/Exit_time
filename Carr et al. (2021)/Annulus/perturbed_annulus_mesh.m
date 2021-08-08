function gmsh = perturbed_annulus_mesh(theta, R_inner, R_outer, gmsh_path, ref, Nbpts)
%AnnulusMeshGeneration(theta, R_inner, R_outer, gmsh_path, ref, Nbpts):
%Constructs the geometry and the mesh for an annulus.
%INPUTS: theta: The angles used to define the geometry.
%        R_inner: A function handle which defines the inner radius.
%        R_outer: A function handle which defines the outer radius.
%        gmsh_path: Where to find gmsh in the current directory
%        ref: Characteristic length for gmsh.
%        Nbpts: Number of points used to capture boundary.
%OUTPUTS: gmsh: This is a structure which contains:
%               gmsh.elements: The connectivity matrix for the elements
%               (triangles) used in the mesh.
%               gmsh.nodes: The coordinates for the nodes in the mesh.
%               gmsh.no_nodes: How many nodes there are in the mesh.
%               gmsh.no_elements: How many elements there are in the mesh.
%               gmsh.outer_boundary_nodes: The indices in gmsh.nodes for
%                 the nodes which are on the outer boundary of the annulus.
%               gmsh.inner_boundary_nodes: The indices in gmsh.nodes for
%                 the nodes which are on the inner boundary of the annulus.
%% Construct Geometry
% Locate File
fid = fopen('Annulus.geo', 'w'); % Open and 'w'rite

% Define parameters
fprintf(fid, 'r = %g;\n', ref); % characteristic length
fprintf(fid,'Mesh.Algorithm = %i;\n', 6);
fprintf(fid,'Mesh.Format = 50;\n'); % matlab output (creates Annulus.m)
% Create Cartesian points
r_outer = R_outer(theta); % outer radii values
r_inner = R_inner(theta); % inner radii values
[x_outer, y_outer] = pol2cart(theta, r_outer);
[x_inner, y_inner] = pol2cart(theta, r_inner);

% Outer boundary points
for k = 1:Nbpts
    fprintf(fid, 'Point(%i) = {%2.15f, %2.15f, 0, r};\n', k, x_outer(k), y_outer(k));
end

% Inner boundary points
for k = (Nbpts+1):2*Nbpts
    fprintf(fid, 'Point(%i) = {%2.15f, %2.15f, 0, r};\n', k, x_inner(k-Nbpts), y_inner(k-Nbpts));
end

% Spline through outer boundary
fprintf(fid,'BSpline(1) = {');
fprintf(fid,'%i,',1:Nbpts);
fprintf(fid,'%i,1};\n',Nbpts);

% Spline through inner boundary
fprintf(fid,'BSpline(2) = {');
fprintf(fid,'%i,',(Nbpts+1):2*Nbpts);
fprintf(fid,'%i,%i};\n',2*Nbpts,Nbpts+1);

% Define plane surface
fprintf(fid,'Line Loop(3) = {1};\n');
fprintf(fid,'Line Loop(4) = {2};\n');
fprintf(fid,'Plane Surface(5) = {3, 4};\n');

% Define physical properties
fprintf(fid,'Physical Line(55) = {1};\n'); % this is the outer line
fprintf(fid,'Physical Line(75) = {2};\n'); % this is the inner line
fprintf(fid,'Physical Surface(99) = {5};\n');

% Close
fclose(fid);

% Generate Mesh
system([gmsh_path, ' Annulus.geo -2']);