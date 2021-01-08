function form = perturbed_ellipse_laplace(a, b, f, N)
%perturbed_ellipse_laplace: Solves del^2 T = 0 on the ellipse x^2/a^2 + y^2/b^2 <= 1,
%where T = f(theta) on x^2/a^2 + y^2/b^2 = 1, where theta is defined by the
%parameter in x = a cos theta, y = b sin theta, 0 <= theta <= 2pi. N terms
%are taken in the Fourier series.
if ~strcmpi(class(a), 'sym') | ~strcmpi(class(b), 'sym')
    if a <= b
        error('a must be greater than b.');
    end
end
syms theta n x y jj kk
An = sym(zeros(1, N+1));
Bn = sym(zeros(1, N+1));
uk = sym(zeros(1, N+1));
vk = sym(zeros(1, N+1));
DU = sym(zeros(1, N+1));
DV = sym(zeros(1, N+1));
Sn = sym(zeros(1, N+1));
Hn = sym(zeros(1, N+1));
Un = sym(zeros(1, N+1));
Vn = sym(zeros(1, N+1));

%% Compute Fourier Coefficients
A0 = 1/sym(pi)*int(f(theta), 0, 2*sym(pi)); % [Eq. 26]
for n = 1:N 
    An(n) = 1/(sym(pi))*int(f(theta)*cos(n*theta), 0, 2*sym(pi)); % [Eq. 26]
    Bn(n) = 1/(sym(pi))*int(f(theta)*sin(n*theta), 0, 2*sym(pi)); % [Eq. 26]
end

%% Compute Re z^k and Im z^k
for k = 1:N
    uk(k) = symsum(nchoosek(k,2*jj)*x^(k-2*jj)*(-1)^jj*y^(2*jj), jj, 0, floor(k/2)); % [Eq. 27]
    vk(k) = symsum(nchoosek(k,2*jj+1)*x^(k-2*jj-1)*(-1)^jj*y^(2*jj+1), jj, 0, floor((k+1)/2)); % [Eq. 28]
end
uk = [1 uk];
vk = [0 vk];
%% Compute DU
for n = 1:2:N
    C = 2*CosinePolynomialCoeff(n);
    for r=0:((n-1)/2)
        DU(n) = DU(n) + C(2*r+2)*(a^2-b^2)^((n-1)/2 - r)*uk(2*r+2); % [Eq. 31]
    end
end
for n =2:2:N
    C = 2*CosinePolynomialCoeff(n);
    for r=0:(n/2)
        DU(n) = DU(n) + C(2*r+1)*(a^2-b^2)^(n/2-r)*uk(2*r+1); % [Eq. 29]
    end
end

%% Compute DV
for n = 1:2:N
    C = 2*CosinePolynomialCoeff(n);
    for r=0:((n-1)/2)
        DV(n) = DV(n) + C(2*r+2)*(a^2-b^2)^((n-1)/2 - r)*vk(2*r+2); % [Eq. 32]
    end
end
for n =2:2:N
    C = 2*CosinePolynomialCoeff(n);
    for r=0:(n/2)
        DV(n) = DV(n) + C(2*r+1)*(a^2-b^2)^(n/2-r)*vk(2*r+1); % [Eq. 30]
    end
end

%% Compute Denominators
for n = 1:N
    Sn(n) = (a+b)^n - (a-b)^n; % [Eq. 29, 31]
    Hn(n) = (a+b)^n + (a-b)^n; % [Eq. 30, 32]
end

%% Compute Un and Vn
for n=1:N
    Un(n) = DU(n) / Hn(n); % [Eq. 29, 31]
    Vn(n) = DV(n) / Sn(n); % [Eq. 30, 32]
end

%% Compute Solution
form = 1/2 * A0; % [Eq. 33]
for n = 1:N
    form = form + An(n)*Un(n)+Bn(n)*Vn(n); % [Eq. 33]
end
function [coeff, poly] = CosinePolynomialCoeff(n)
%CosinePolynomialCoeff(n): Outputs the coefficients in the expansion of cos(nx)
%in terms of cos(x) in ascending order, e.g. cos(7x) = -7 cos(x) + 56
%cos^3(x) - 112 cos^5(x) + 64 cos^7(x), when coeff = [0 -7 0 56 0 -112 0 64].
% poly is the actual expansion, where c stands for cos(x).
syms c
coeff = fliplr(coeffs(chebyshevT(n, c), 'all'));
poly = poly2sym(fliplr(coeff), c);