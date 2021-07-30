function ExitTime = unperturbed_ellipse_walk(x, y, Delta, tau, P, N, a, b)
%unperturbed_ellipse_walk: Simulates a random walk on an ellipse N times and computes the
%mean exit time.
%INPUTS: [x, y]: The starting positions.
%        Delta: The step-size.
%        tau: The step duration.
%        P: The probability that the walk moves in each step duration.
%        N: The number of simulations.
%        a: Width of ellipse.
%        b: Height of ellipse.
%OUTPUT: ExitTime: The mean exit time for the walk.
% Initiate
ExitTime = zeros(length(x), 1);
for n=1:N
    for i = 1:size(x, 1)
        k = 1;
        while x(i, k)^2/a^2 + y(i, k)^2/b^2 < 1 
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