function T = perturbed_disc_perturbation(R, D, n, N, g, epsilon)
%perturbed_disc_perturbation: Finds the perturbation solution T(r, theta) = T_0 +
%epsilon T_1 + epsilon^2 T_2 + ... + epsilon^n T_n to the problem of first
%passage time in a perturbed geometry with boundary {R^*(theta) : 0 <=
%theta < 2pi}, where R^*(theta) = R(1+ep g(theta)), where epsilon is a
%small parameter, R is the radius of the unperturbed boundary, and g(theta)
%is a smooth 2pi periodic function.
%The diffusivity is given as D, and we are taking N terms in the Fourier
%series with a parameter epsilon << 1.
syms r; assume(r > 0);
syms theta;
syms k; assume(k, 'integer'); assumeAlso(k > 0);
syms T0 % The reason to avoid using symbolic FUNCTIONS is for ease of use in the loop.
syms dT0
syms T [1 n] % main function
T = [T0 T];
T(1) = (R^2 - r^2)/(4*D); % [Equation (4)]
syms Hj
for j=1:n
    Hj = 0;
    for p=1:j
        Hj = Hj - (R*g(theta))^p/factorial(p)*subs(diff(T(j-p+1), r, p), r, R); % [Equation (13)]
    end
    A0 = 1/sym(pi) * intfix(Hj, theta, 0, 2*sym(pi)); % The reason to use a different function is because MATLAB often fails to integrate orthogonal products symbolically. [Equation (15)]
    Ak = 1/(sym(pi)*R^k) * intfix(Hj*cos(k*theta), theta, 0, 2*sym(pi)); % [Equation (15)]
    Bk = 1/(sym(pi)*R^k) * intfix(Hj*sin(k*theta), theta, 0, 2*sym(pi)); % [Equation (15)]
    [~, dAk] = numden(Ak);
    [~, dBk] = numden(Bk);
    d1 = solve(dAk == 0, k);
    d2 = solve(dBk == 0, k);
    divideZeroIdx = double(union(d1, d2))';
    divideZeroIdx = divideZeroIdx(divideZeroIdx > 0 & divideZeroIdx <= N);
    T(j+1) = 1/2 * A0; % [Equation (14)]
    if ~isempty(divideZeroIdx)
        for kk = divideZeroIdx
            Ak1 = 1/(sym(pi)*R^kk) * intfix(Hj*cos(kk*theta), theta, 0, 2*sym(pi)); % [Equation (15)]
            Bk1 = 1/(sym(pi)*R^kk) * intfix(Hj*sin(kk*theta), theta, 0, 2*sym(pi)); % [Equation (15)]
            T(j+1) = T(j+1) + subs(r^k*(Ak1*cos(k*theta)+Bk1*sin(k*theta)), k, kk); % [Equation (15)]
        end
    end
    for kk = setdiff(1:N, divideZeroIdx)
        T(j+1) = T(j+1) + subs(r^k*(Ak*cos(k*theta)+Bk*sin(k*theta)), k, kk); % [Equation (14)]
    end
end
T = sum(epsilon.^(0:n) .* T); % [Equation (9)]
function I = intfix(expr, var, a, b)
I = int(expr, var);
I = simplifyFraction(subs(I, b) - subs(I, a));