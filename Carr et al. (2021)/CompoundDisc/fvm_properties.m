function mesh = fvm_properties(gmsh, R_vec)
% Computes mesh properties (normal vectors, control volume areas, gradient
% approximation coefficients etc)
%
% Inputs:
% nodes (matrix of size: no_nodes by 2)
%     Matrix containing the coordinates of the nodes of the mesh such that nodes(i,1) and nodes(i,2)
%     denote the the x and y coordinates of the ith node.
% elements (matrix of size: no_elements by 3)
%     Matrix containing the vertices of each element, such that elements(k,:) returns the three 
%     vertex/nodes of the kth element stored counter-clockwise. For example, if the 3rd element has 
%     vertices corresponding to nodes 12, 73 and 106 (order couter-clockwise) then elements(3,:) 
%     returns the vector [12, 73, 106].
%
% Outputs:
% mesh (structure array)
%     Structure array containing the geometrical properties of the mesh (normal vectors, control 
%     volume areas, gradient approximation coefficients etc).
%
% Created by: Dr Elliot Carr (2019), Queensland University of Technology.
% Modified by: Josh Wilson (2020), Queensland University of Technology.

% Get the number of nodes and control volumes from the mesh object
no_elements = size(gmsh.elements, 1);
no_nodes = size(gmsh.nodes, 1);

% Create space to store node information
variable = zeros(no_nodes,1);
mapping = zeros(no_nodes,1);

% Solution on boundary known from Dirichlet boundary conditions
no_variables = 0;
for i = 1:no_nodes
    if ismember(i, gmsh.boundary_nodes) 
        % Known from Dirichlet boundary condition
        variable(i) = 0;
    else
        % Variables are nodes not situated on the boundary 
        variable(i) = 1;
    end
    mapping(i) = i;
end

% Initialise
V  = zeros(no_nodes,1);
nx = zeros(no_elements,3);
ny = zeros(no_elements,3);
s  = zeros(no_elements,6);
E  = zeros(no_elements,1);
L  = zeros(no_elements,3);
layer = zeros(no_elements,1);
cell_area = 0;

% Compute geometrical properties of mesh
for k = 1:no_elements
    
    v = gmsh.elements(k,:); % vertices of kth element
    xv   = gmsh.nodes(v,1); % x-coordinates of the three nodes/vertices
    yv   = gmsh.nodes(v,2); % y-coordinates of the three nodes/vertices
    
    % Element centroid
    centroid = [sum(xv) / 3, sum(yv) / 3];
    
    % Element edge midpoint
    m1   = [(xv(1) + xv(2)) / 2, (yv(1) + yv(2)) / 2];
    m2   = [(xv(2) + xv(3)) / 2, (yv(2) + yv(3)) / 2];
    m3   = [(xv(1) + xv(3)) / 2, (yv(1) + yv(3)) / 2];
    
    % Control volume edge vectors
    e1 = centroid - m1;
    e2 = centroid - m2;
    e3 = centroid - m3;
    
    % Control volume edge lengths
    L(k,1) = norm(e1,2);
    L(k,2) = norm(e2,2);
    L(k,3) = norm(e3,2);
    
    % Unit normals
    nx(k,1) = e1(2)/L(k,1); ny(k,1) = -e1(1)/L(k,1);
    nx(k,2) = e2(2)/L(k,2); ny(k,2) = -e2(1)/L(k,2);
    nx(k,3) = e3(2)/L(k,3); ny(k,3) = -e3(1)/L(k,3);
    
    % Vectors connecting nodes to centroid
    p1 = centroid - [xv(1), yv(1)];
    p2 = centroid - [xv(2), yv(2)];
    p3 = centroid - [xv(3), yv(3)];
    
    % Vectors connecting midpoints
    q1 = m1 - m3;
    q2 = m2 - m1;
    q3 = m3 - m2;
    
    % Sub-control volume areas
    S1 = 1/2 * norm(cross([p1, 0], [q1, 0]));
    S2 = 1/2 * norm(cross([p2, 0], [q2, 0]));
    S3 = 1/2 * norm(cross([p3, 0], [q3, 0]));
    
    % Control volume areas
    V(v(1)) = V(v(1)) + S1;
    V(v(2)) = V(v(2)) + S2;
    V(v(3)) = V(v(3)) + S3;
     
    % Element areas
    E(k) = S1 + S2 + S3;
    
    % Unit cell area
    cell_area = cell_area + E(k);
    
    % Shape function information
    detA = xv(2)*yv(3) - yv(2)*xv(3) - xv(1)*yv(3) + xv(3)*yv(1) + xv(1)*yv(2) - xv(2)*yv(1);
    s(k,1) = (yv(2) - yv(3)) / detA;
    s(k,2) = (yv(3) - yv(1)) / detA;
    s(k,3) = (yv(1) - yv(2)) / detA;
    s(k,4) = (xv(3) - xv(2)) / detA;
    s(k,5) = (xv(1) - xv(3)) / detA;
    s(k,6) = (xv(2) - xv(1)) / detA;
    
    [theta, r] = cart2pol(centroid(1), centroid(2));
    for i = 1:length(R_vec) - 1
        if R_vec{i}(theta) < r && r < R_vec{i+1}(theta)
            layer(k) = i;
            break
        end
    end
    
end

mesh.elements     = gmsh.elements;
mesh.nodes        = gmsh.nodes;
mesh.no_elements  = no_elements;
mesh.no_nodes     = no_nodes;
mesh.no_variables = no_variables;
mesh.nx           = nx;
mesh.ny           = ny;
mesh.s            = s;
mesh.mapping      = mapping;
mesh.variable     = variable;
mesh.V            = V;
mesh.cell_area    = cell_area;
mesh.E            = E;
mesh.L            = L;
mesh.layer        = layer;

end