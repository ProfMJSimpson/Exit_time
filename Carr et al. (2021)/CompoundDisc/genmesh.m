function gmsh = genmesh(R_vec, meshrefine, numbound, gmsh_path)

% Starting positions
theta = linspace(0, 2*pi, numbound);
gmsh = starting_positions(theta, R_vec, gmsh_path, meshrefine);

