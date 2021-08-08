function gmsh = starting_positions(theta,R_vec,gmsh_path,ref)

Nbpts = length(theta);
Nlays = length(R_vec) - 2;

fid = fopen('mesh.geo', 'w');
fprintf(fid,'r = %g;\n',ref);
fprintf(fid,'Mesh.Algorithm = 6;\n');
fprintf(fid,'Mesh.Format = 50;\n');

r = zeros(Nlays, Nbpts);
x = zeros(Nlays, Nbpts);
y = zeros(Nlays, Nbpts);
count = 0;
for i = 1:Nlays
    r(i, :) = R_vec{i + 1}(theta);
    [x(i, :), y(i, :)] = pol2cart(theta, r(i, :));
    for k = 1:Nbpts
        fprintf(fid,'Point(%i) = {%g,%g,0,r};\n', count + k, x(i, k), y(i, k));
    end
    fprintf(fid,'BSpline(%i) = {',i);
    fprintf(fid,'%i,', count + 1:count + Nbpts);
    fprintf(fid,'%i,%i};\n', count + Nbpts, count + 1);
    fprintf(fid,'Line Loop(%i) = {%i};\n', i, i);
    fprintf(fid,'Physical Line(%i) = {%i};\n', 98 + i, i);
    if(i ~= 1)
        fprintf(fid,'Plane Surface(%i) = {%i, %i};\n', i, i, i - 1);
    else
        fprintf(fid,'Plane Surface(1) = {1};\n');
    end
    count = count + Nbpts;
end

fprintf(fid,'Physical Surface(1) = {');
if(Nlays > 1)
    fprintf(fid,'%i,', 1:Nlays - 1);
end
fprintf(fid,'%i};\n', Nlays);
fclose(fid);
system([gmsh_path,' mesh.geo -2']);

%% Read gmsh .msh file
mesh;
gmsh.nodes = msh.POS(:,1:2);
gmsh.elements = msh.TRIANGLES(:,1:3);
if(Nlays == 1)
    gmsh.boundary_nodes = msh.LINES(msh.LINES(:,3) == 99, 1);
    gmsh.interface_nodes = msh.LINES(msh.LINES(:,3) ~= 99, 1);
else
    gmsh.boundary_nodes = msh.LINES(msh.LINES(:,3) ~= 99, 1);
    gmsh.interface_nodes = msh.LINES(msh.LINES(:,3) == 99, 1);
end
