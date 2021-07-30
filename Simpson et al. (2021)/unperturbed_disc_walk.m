function ExitTime = unperturbed_disc_walk(x, y, R, Delta, tau, P, N, nodes, boundary_nodes)
%unperturbed_disc_walk: Simulates a random walk on a disc for walks
%starting at positions designated by [x, y] until the walk leaves the disc.
%Delta is the step-size, tau is the time-step, P is the probability that
%the walk takes a step during each time epoch, and R is the radius of the
%disc. The points [x y] are in Cartesian coordinates. N is the number of
%realisiations used to compute the mean exit time. 
arguments
    x (:, 1) = 0
    y (:, 1) = 0 % specify as column vectors
    R (1, 1) = 1
    Delta (1, 1) = 1
    tau (1, 1) = 1
    P (1, 1) = 1
    N (1, 1) = 50
    nodes = [x y]
    boundary_nodes = NaN
end
%% Initiate
num_walks = length(x);
ExitTime = zeros(num_walks, N);
walks_inside = PointInside(x, y, R);
walks_inside(boundary_nodes) = false;

%% Main loop
k = 1;
for n=1:N
    while sum(walks_inside) > 0
        move_prob = rand(num_walks, 1);
        theta = 2*pi*rand(num_walks, 1);
        move_walk = move_prob <= P;
        move_idx_1 = walks_inside(1:end);
        move_idx_2 = move_idx_1 & move_walk;
        x(move_idx_2) = x(move_idx_2) + Delta * cos(theta(move_idx_2));
        y(move_idx_2) = y(move_idx_2) + Delta * sin(theta(move_idx_2));
        ExitTime(move_idx_1, k) = ExitTime(move_idx_1, k) + tau;
        walks_inside(move_idx_1)  = PointInside(x(move_idx_1), y(move_idx_1), R);
    end
    x = nodes(:, 1); y = nodes(:, 2);
    walks_inside = PointInside(x, y, R);
    k = k + 1;
end
ExitTime = mean(ExitTime, 2);

% Test if inside boundary
function InsideShape = PointInside(x, y, R)
[~, r] = cart2pol(x, y);
InsideShape = r < R;