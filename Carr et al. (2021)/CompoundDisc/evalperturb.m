function plotvecs = evalperturb(x0, y0, R1star, STfunc, epsilon)
    
    syms theta r

    [tvals, rvals] = cart2pol(x0, y0);

    STmatlab = cell(size(STfunc));
    for i = 1:size(STfunc, 1)
        STmatlab{i, 1} = matlabFunction(STfunc{i, 1}, 'vars', [r, theta]);
        STmatlab{i, 2} = matlabFunction(STfunc{i, 2}, 'vars', [r, theta]);
    end

    pertsol = zeros(length(x0), size(STfunc, 1));
    for j = 1:size(pertsol, 2)
        Svals = STmatlab{j, 1}(rvals, tvals);
        Tvals = STmatlab{j, 2}(rvals, tvals);
        for i = 1:size(pertsol, 1)
            if rvals(i) < R1star(tvals(i))
                pertsol(i, j) = Svals(i);
            else
                pertsol(i, j) = Tvals(i);
            end
        end
    end

    plotvecs = zeros(size(pertsol));
    for i = 0:size(pertsol, 2) - 1
        for j = i + 1:size(pertsol, 2)
            plotvecs(:, j) = plotvecs(:, j) + epsilon^i * pertsol(:, i + 1);
        end
    end

end