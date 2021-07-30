function unperturbed_ellipse_mesh(x, y, gmsh_path, ref, Nbpts)
%unperturbed_ellipse_mesh: Generates the mesh on an unperturbed ellipse
%with boundary values specified by x and y, using the directory gmsh_path
%to locate gmsh. ref is the refinement parameter, and Nbpts is the number
%of points used to capture the boundary.

fid = fopen('Ellipse.geo', 'w');
fprintf(fid,'r = %g;\n',ref);
fprintf(fid,'Mesh.Algorithm = %i;\n', 6);
fprintf(fid,'Mesh.Format = 50;\n'); % matlab output (creates Ellipse.m)
for k = 1:Nbpts
    fprintf(fid,'Point(%i) = {%2.15f,%2.15f,0,r};\n',k,x(k),y(k));
end
fprintf(fid,'BSpline(1) = {');
fprintf(fid,'%i,',1:Nbpts);
fprintf(fid,'%i,1};\n',Nbpts);
fprintf(fid,'Line Loop(1) = {1};\n');
fprintf(fid,'Physical Line(99) = {1};\n');
fprintf(fid,'Plane Surface(1) = {1};\n');
fprintf(fid,'Physical Surface(1) = {1};\n');
fclose(fid);
system([gmsh_path,'Ellipse.geo -2']);