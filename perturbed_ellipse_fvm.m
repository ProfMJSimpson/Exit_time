function T = perturbed_ellipse_fvm(D, nodes, elements, boundary_nodes)
%perturbed_ellipse_fvm: Solves the PDE del^2 T = -1/D on an ellipse.
%INPUTS: D: The diffusivity.
%        nodes: The nodes to solve at.
%        elements: The connectivity matrix for the triangular mesh.
%        boundary_nodes: The boundary nodes on the ellipse.
%OUTPUS: T: The solution to the PDE at each node.r
no_elements = size(elements, 1); % number of elements
no_nodes = size(nodes, 1); % number of nodes


% Specify nodes on boundary
variable = zeros(no_nodes, 1);
for i = 1:no_nodes
    if ismember(i, boundary_nodes)
        variable(i) = 0;
    else
        variable(i) = 1;
    end
end

% Initialise
V  = zeros(no_nodes,1); % volumes
nx = zeros(no_elements,3); % x-coordinate normal
ny = zeros(no_elements,3); % y-coordinate normal
s  = zeros(no_elements,6); % for the gradient approximation
E  = zeros(no_elements,1); % area of each element
L  = zeros(no_elements,3); % length of each segment
cell_area = 0; % area of cell

% Compute geometrical properties of mesh
for k = 1:no_elements
    
    v = elements(k,:); % vertices of kth element
    xv   = nodes(v,1); % x-coordinates of the three nodes/vertices
    yv   = nodes(v,2); % y-coordinates of the three nodes/vertices
    
    % Element centroid
    centroid = [mean(xv), mean(yv)];
    
    % Element edge midpoint
    m1   = [mean(xv([1 2])), mean(yv([1 2]))];
    m2   = [mean(xv([2 3])), mean(yv([2 3]))];
    m3   = [mean(xv([1 3])), mean(yv([1 3]))];
    
    % Control volume edge vectors
    e1 = centroid-m1;
    e2 = centroid-m2;
    e3 = centroid-m3;
    
    % Control volume edge lengths
    L(k,1) = norm(e1,2);
    L(k,2) = norm(e2,2);
    L(k,3) = norm(e3,2);
    
    % Unit normals
    nx(k,1) = e1(2)/L(k,1); ny(k,1) = -e1(1)/L(k,1); % To rotate a vector 90 degrees **CW**, (x, y) -> (y, -x)
    nx(k,2) = e2(2)/L(k,2); ny(k,2) = -e2(1)/L(k,2);
    nx(k,3) = e3(2)/L(k,3); ny(k,3) = -e3(1)/L(k,3);
    
    % Vectors connecting nodes to centroid
    p1 = centroid-[xv(1) yv(1)];
    p2 = centroid-[xv(2) yv(2)];
    p3 = centroid-[xv(3) yv(3)];
    
    % Vectors connecting midpoints
    q1 = m1-m3;
    q2 = m2-m1;
    q3 = m3-m2;
    
    % Sub-control volume areas
    S1 = norm(p1, 2)*norm(q1, 2)/2;
    S2 = norm(p2, 2)*norm(q2, 2)/2;
    S3 = norm(p3, 2)*norm(q3, 2)/2;
    
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
    s(k,1) = (yv(2)-yv(3)) / detA;
    s(k,2) = (yv(3)-yv(1)) / detA;
    s(k,3) = (yv(1)-yv(2)) / detA;
    s(k,4) = (xv(3)-xv(2)) / detA;
    s(k,5) = (xv(1)-xv(3)) / detA;
    s(k,6) = (xv(2)-xv(1)) / detA;
    
end

% Initialise
A = sparse(no_nodes, no_nodes);
b = -V/D; % [Eq. 1]

% Loop over elements and build up A
for k = 1:no_elements
    v = elements(k, :); % vertices of element k
    for j = 1:3 % loop over 3 edges in element k
        if j == 1 % control volume edge j (shared by vertex j and jnb)
            jnb = 2;
        elseif j == 2
            jnb = 3;
        elseif j == 3
            jnb = 1;
        end
        
        % Add contribution to node at vertex j
        A(v(j),v(1)) = A(v(j),v(1)) + (s(k,1)*nx(k,j) + s(k,4)*ny(k,j))*L(k,j);
        A(v(j),v(2)) = A(v(j),v(2)) + (s(k,2)*nx(k,j) + s(k,5)*ny(k,j))*L(k,j);
        A(v(j),v(3)) = A(v(j),v(3)) + (s(k,3)*nx(k,j) + s(k,6)*ny(k,j))*L(k,j);
        
        % Add contribution to node at vertex jnb
        A(v(jnb),v(1)) = A(v(jnb),v(1)) - (s(k,1)*nx(k,j) + s(k,4)*ny(k,j))*L(k,j);
        A(v(jnb),v(2)) = A(v(jnb),v(2)) - (s(k,2)*nx(k,j) + s(k,5)*ny(k,j))*L(k,j);
        A(v(jnb),v(3)) = A(v(jnb),v(3)) - (s(k,3)*nx(k,j) + s(k,6)*ny(k,j))*L(k,j);
    end
end

% Discrete equations for boundary

for i = 1:no_nodes
    if variable(i) == 0 % [Eq. 2]
        A(i, i) = 1;
        A(i, [1:i-1 i+1:no_nodes]) = 0;
        b(i) = 0;
    end
end

% Solve linear system
T = A \ b;

