function T = perturbed_ellipse_perturbation(a, b, g, h, D, perturb_degree, N, epsilon)
%perturbed_ellipse_perturbation: Finds the perturbation solution T = T0 + epsilon T1 +
%epsilon^2 T2 + ... + epsilon^n Tn to the problem of mean first passage
%time in a perturbed ellipse with boundary parametrised as x = a(1 +
%epsilon g(theta))cos(theta), y = a(1 + epsilon h(theta))sin(theta), 0 <=
%theta <= 2pi. Here, a > b.
syms x y
syms theta
syms n; assume(n, 'integer'); assumeAlso(n > 0);
m = perturb_degree + 1;
syms T [1 m]
syms H(theta)
temp_err = warning('error', 'symbolic:solve:warnmsg3');
if ~strcmpi(class(a), 'sym') | ~strcmpi(class(b), 'sym')
    if a < b
        error('a must be greater than b.')
    end
end
T(1) = a^2*b^2/(2*D*(a^2+b^2))*(1-x^2/a^2-y^2/b^2); % [Equation 6]
for ell = 1:perturb_degree
    H(theta) = 0;
    for k=1:ell
        for i=0:k
            H(theta) = H(theta) - 1/factorial(k) * nchoosek(k, i) * subs(diff(diff(T(ell-k+1), x, i), y, k-i), [x y], [a*cos(theta) b*sin(theta)]) * a^i*b^(k-i)*g(theta).^i.*h(theta)^(k-i).*cos(theta)^i.*sin(theta).^(k-i); % [Equation 20]
        end
    end
    T(ell+1) = perturbed_ellipse_laplace(a, b, H, N); % [Eq. 24]
end
T = dot(epsilon.^(0:perturb_degree), T); % [Eq. 18]
warning(temp_err);