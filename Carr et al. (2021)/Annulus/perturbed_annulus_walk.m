%function ExitTime = perturbed_annulus_walk(x, y, Delta, tau, P, N, nodes, outer_boundary_nodes, inner_boundary_nodes, inner_type)
function ExitTime = perturbed_annulus_walk(x, y, Delta, tau, P, N, R_inner, R_outer, inner_type)
%ExitTime: Simulates a random walk on an annulus N times and computes the
%mean exit time.
%INPUTS: [x, y]: The starting positions.
%        Delta: The step-size.
%        tau: The step duration.
%        P: The probability that the walk moves in each step duration.
%        N: The number of simulations.
%        nodes: The nodes used in the mesh which define [x, y].
%        outer_boundary_nodes: The indices in nodes which give the outer
%               boundary of the annulus.
%        inner_boundary_nodes: The indices in nodes which give the inner
%               boundary of the annulus.
%        inner_type: 'reflect' if the inner boundary should be reflective,
%               and 'absorb' if the inner boundary should absorb.
%OUTPUT: ExitTime: The mean exit time for the walk.
%NOTES: By default, the outer boundary absorbs.

% Initiate
num_walks = length(x);
ExitTime = zeros(num_walks, 1);

% Main loop
for n = 1:N % realisation loop
    
    % Walk loop
    for i = 1:num_walks % walk loop
        k = 1;
        
        % Perform walk according to reflective inner boundary
        if strcmpi(inner_type, 'reflect')
            
            % Perform walk until it leaves outer boundary
            %while PointInside(x(i, k), y(i, k), nodes(outer_boundary_nodes, 1), nodes(outer_boundary_nodes, 2))
            while PointInside(x(i, k), y(i, k), R_outer)
                move_prob = rand;
                if move_prob <= P
                    theta = 2*pi*rand; % simulate random angle
                    x(i, k+1) = x(i, k) + Delta * cos(theta);
                    y(i, k+1) = y(i, k) + Delta * sin(theta);
                    %if PointInside(x(i, k+1), y(i, k+1), nodes(inner_boundary_nodes, 1), nodes(inner_boundary_nodes, 2))
                    if PointInside(x(i, k+1), y(i, k+1), R_inner)
                        x(i, k+1) = x(i, k);
                        y(i, k+1) = y(i, k); % reflect back
                    end
                else
                    x(i, k+1) = x(i, k);
                    y(i, k+1) = y(i, k);
                end
                ExitTime(i) = ExitTime(i) + tau;
                k = k + 1;
            end
            
        elseif strcmpi(inner_type, 'absorb')
            
            % Perform walk until it leaves either boundary
            %while PointInside(x(i, k), y(i, k), boundary_nodes_x, boundary_nodes_y)
            while PointInside(x(i, k), y(i, k), R_outer) && ~PointInside(x(i, k), y(i, k), R_inner)
                move_prob = rand;
                if move_prob <= P
                    theta = 2*pi*rand;
                    x(i, k+1) = x(i, k) + Delta * cos(theta);
                    y(i, k+1) = y(i, k) + Delta * sin(theta);
                else
                    x(i, k+1) = x(i, k);
                    y(i, k+1) = y(i, k);
                end
                ExitTime(i) = ExitTime(i) + tau;
                k = k + 1;
            end
        end
    end
    
    clc
    disp(['Random Walk completion: ',num2str(round(100*n/N)),'%']);
end
ExitTime = ExitTime/N;

% Test if inside boundary
% function InsideShape = PointInside(x, y, boundary_x, boundary_y)
% InsideShape = inpolygon(x, y, boundary_x, boundary_y);
% end
% end
function InsideShape = PointInside(x, y, R)
[theta, r] = cart2pol(x, y);
InsideShape = r < R(theta);
