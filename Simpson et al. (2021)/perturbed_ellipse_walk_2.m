function ExitTime = perturbed_ellipse_walk_2(x, y, Delta, tau, P, N, ...
    nodes,boundary_nodes)
%perturbed_ellipse_walk_2: Simulates a random walk on an ellipse N times and computes the
%mean exit time.
%INPUTS: [x, y]: The starting positions.
%        Delta: The step-size.
%        tau: The step duration.
%        P: The probability that the walk moves in each step duration.
%        N: The number of simulations.
%        x_fnc: x-perturbation function.
%        y_fnc: y-perturbation function.
%OUTPUT: ExitTime: The mean exit time for the walk.
%NOTES: If you encounter errors in inpoly2, you must add 
% if isempty(ddxy)
%        return
%    end
% to the part before the statement if 
% if (ddxy(1) > ddxy(2))
%
% Initiate
num_walks = length(x);
ExitTime = zeros(num_walks, 1);
boundary_nodes_x = nodes(boundary_nodes, 1)';
boundary_nodes_y = nodes(boundary_nodes, 2)';
walks_inside = PointInside(x, y, boundary_nodes_x, boundary_nodes_y);
% Main loop
for n=1:N
    while sum(walks_inside) > 0
        move_prob = rand(num_walks, 1);
        theta = 2*pi*rand(num_walks, 1);
        move_walk = move_prob <= P;
        move_idx_1 = walks_inside(1:end);
        move_idx_2 = move_idx_1 & move_walk;
        x(move_idx_2) = x(move_idx_2) + Delta * cos(theta(move_idx_2));
        y(move_idx_2) = y(move_idx_2) + Delta * sin(theta(move_idx_2));
        ExitTime(move_idx_1) = ExitTime(move_idx_1) + tau;
        walks_inside(move_idx_1)  = PointInside(x(move_idx_1), y(move_idx_1), boundary_nodes_x, boundary_nodes_y);
    end
    x = nodes(:, 1); y = nodes(:, 2);
    walks_inside = PointInside(x, y, boundary_nodes_x, boundary_nodes_y);
end
ExitTime = ExitTime/N;

% Test if inside boundary
function InsideShape = PointInside(x, y, boundary_x, boundary_y)
if x <= min(boundary_x') | x >= max(boundary_x') | y <= min(boundary_y') | y >= max(boundary_y')
    InsideShape = false;
else
    InsideShape = reshape(inpoly2([x(:) y(:)], [boundary_x' boundary_y']), size(x));
end