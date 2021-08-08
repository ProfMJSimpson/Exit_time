function u = fvm_heterogeneous(D, mesh)
% Solves the homogenization boundary value problems
%
% Inputs:
% D (double)
%     Global diffusivity for the problem
% mesh (structure array)
%     Structure array containing the geometrical properties of the mesh (normal vectors, control 
%     volume areas, gradient approximation coefficients etc).
%
% Outputs:
% u1 (vector of size: no_nodes by 1)
%     Finite volume solution of homogenization boundary value problem required for Deff(:,1)
% u2 (vector of size: no_nodes by 1)
%     Finite volume solution of homogenization boundary value problem required for Deff(:,2)
%
% Dr Elliot Carr (2019), Queensland University of Technology.

elements        = mesh.elements;
no_elements     = mesh.no_elements;
no_nodes        = mesh.no_nodes;
nx              = mesh.nx;
ny              = mesh.ny;
s               = mesh.s;
mapping         = mesh.mapping;
variable        = mesh.variable;
V               = mesh.V;
L               = mesh.L;
layer           = mesh.layer;

% Initialise
A = sparse(no_nodes,no_nodes);
b = sparse(no_nodes,1);

% Loop over elements and build up A and b
for k = 1:no_elements
    
    v = elements(k,:); % vertices of element k
    vn = mapping(v); % mapping to account for periodic conditions
    
    for j = 1:3 % loop over 3 edges in element k

        if j == 1 % control volume edge j (shared by vertex j and jnb)
            jnb = 2;
        elseif j == 2
            jnb = 3;
        elseif j == 3
            jnb = 1;
        end
        
        % Add contribution to node at vertex j
        A(vn(j),vn(1)) = A(vn(j),vn(1)) + D(layer(k))*(s(k,1)*nx(k,j) + s(k,4)*ny(k,j))*L(k,j);
        A(vn(j),vn(2)) = A(vn(j),vn(2)) + D(layer(k))*(s(k,2)*nx(k,j) + s(k,5)*ny(k,j))*L(k,j);
        A(vn(j),vn(3)) = A(vn(j),vn(3)) + D(layer(k))*(s(k,3)*nx(k,j) + s(k,6)*ny(k,j))*L(k,j);
        
        % Add contribution to node at vertex jnb
        A(vn(jnb),vn(1)) = A(vn(jnb),vn(1)) - D(layer(k))*(s(k,1)*nx(k,j) + s(k,4)*ny(k,j))*L(k,j);
        A(vn(jnb),vn(2)) = A(vn(jnb),vn(2)) - D(layer(k))*(s(k,2)*nx(k,j) + s(k,5)*ny(k,j))*L(k,j);
        A(vn(jnb),vn(3)) = A(vn(jnb),vn(3)) - D(layer(k))*(s(k,3)*nx(k,j) + s(k,6)*ny(k,j))*L(k,j);
        
    end
    
end

% Discrete equations for boundary nodes
for i = 1:no_nodes
    if variable(i) == 0 % Solution known from Dirichlet BC
        A(i, :) = zeros(1, no_nodes);
        A(i,i) = 1;
        b(i) = 0;
    else
        b(i) = -V(i);
    end
end

% Solve linear systems
u = A\b; % Solution for first escape