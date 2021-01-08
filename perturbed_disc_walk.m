function ExitTime = perturbed_disc_walk(x, y, Delta, tau, P, N, R)
%perturbed_disc_walk: Simulates a random walk on a disc for walks starting at
%positions designated by [x, y] until the walk leaves the disc. Delta is
%the step-size, tau is the time-step, P is the probability that the walk
%takes a step. The points [x y] are in Cartesian coordinates. 
arguments
    x (:, 1) = 0
    y (:, 1) = 0 % specify as column vectors
    Delta (1, 1) = 1
    tau (1, 1) = 1
    P (1, 1) = 1
    N = 50
    R = NaN;
end
ExitTime = zeros(length(x), 1);
for n = 1:N
    for i = 1:size(x, 1)
        k = 1;
        while PointInside(x(i, k), y(i, k), R)
            rand_num = rand;
            if rand_num <= P
                theta = 2*pi*rand; % simulate random angle
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
ExitTime = ExitTime/N;
function InsideShape = PointInside(x, y, R)
[theta, r] = cart2pol(x, y);
InsideShape = r < R(theta);
