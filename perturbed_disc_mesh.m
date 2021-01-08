function perturbed_disc_mesh(theta,R,gmsh_path,ref,Nbpts)
%perturbed_disc_mesh: Generates the mesh on a perturbed disc with boundary
%function R using the angles specified in theta to capture the boundary.
%The gmsh_path argument indicates where in the working directory MATLAB may
%find gmsh, ref is used to control the number of starting positions, and
%Nbpts is the number of boundary points.

fid = fopen('mesh.geo', 'w');
fprintf(fid,'r = %g;\n',ref);
fprintf(fid,'Mesh.Algorithm = %g;\n',6);
fprintf(fid,'Mesh.Format = 50;\n'); % matlab output (creates mesh.m)
r = R(theta);
[x,y] = pol2cart(theta,r);
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
system([gmsh_path,'mesh.geo -2']);