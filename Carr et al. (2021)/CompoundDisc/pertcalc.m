function STfunc = pertcalc(R1, R2, D, Gfunc, maxorder, maxfourier, STfunc_old)

    syms r theta n
    assume(r >= 0);
    assumeAlso(n >= 0)

    D1 = D(1);
    D2 = D(2);
    
    g = Gfunc(theta);
    gdash = diff(g, theta);
    
    if(nargin < 6)
        error('Not enough input arguments for perturbsol.m')
    elseif(nargin ~= 7)
        STfunc = cell(maxorder + 1, 1);
        STfunc{1, 1} = (R1^2 - r.^2)/(4*D1) + (R2^2 - R1^2)/(4*D2);
        STfunc{1, 2} = (R2^2 - r.^2)/(4*D2);
        initorder = 1;
    else
        STfunc = cell(max(length(STfunc_old), maxorder + 1), 2);
        for i = 1:size(STfunc_old, 1)
            STfunc{i, 1} = STfunc_old{i, 1};
            STfunc{i, 2} = STfunc_old{i, 2};
        end
        initorder = size(STfunc_old, 1);
    end
    
    for order = initorder:maxorder
        f1 = 0;
        for k = 1:order
            f1 = f1 - (R1 * g).^k / factorial(k) ...
                .* (subs(diff(STfunc{order - k + 1, 1}, r, k), r, R1) ...
                - subs(diff(STfunc{order - k + 1, 2}, r, k), r, R1));
        end
        
        f2 = 0;
        for k = 1:order
            f2 = f2 - (R1 * g).^k / factorial(k) ...
                .* (D1*subs(diff(STfunc{order - k + 1, 1}, r, k+1), r, R1) ...
                - D2*subs(diff(STfunc{order - k + 1, 2}, r, k+1), r, R1));
        end
        for k = 1:order - 1
            f2 = f2 - 2*g .* (R1 * g).^k / factorial(k) ...
                .* (D1*subs(diff(STfunc{order - k, 1}, r, k+1), r, R1) ...
                - D2*subs(diff(STfunc{order - k, 2}, r, k+1), r, R1));
        end
        for k = 1:order - 2
            f2 = f2 - g.^2 .* (R1 * g).^k / factorial(k) ...
                .* (D1*subs(diff(STfunc{order - k - 1, 1}, r, k+1), r, R1) ...
                - D2*subs(diff(STfunc{order - k - 1, 2}, r, k+1), r, R1));
        end
        for k = 1:order - 1
            f2 = f2 + gdash/R1 .* (R1 * g).^k / factorial(k) ...
                .* (D1*subs(diff(diff(STfunc{order - k, 1}, r, k), theta), r, R1) ...
                - D2*subs(diff(diff(STfunc{order - k, 2}, r, k), theta), r, R1));
        end
        
        Dstar = n*D1*D2*(R1^(n-1) + R2^(2*n)*R1^(-n-1)) / ((D2 - D1)*R1^n + (D2 + D1)*R2^(2*n)*R1^(-n));

        v0 = D2 / (R1 * log(R1 / R2)) * int(f2, theta, 0, 2*pi);
        v1 = Dstar * (int(f1*cos(n*theta), theta, 0, 2*pi) - R1 / (D1 * n) * int(f2*cos(n*theta), theta, 0, 2*pi));
        v2 = Dstar * (int(f1*sin(n*theta), theta, 0, 2*pi) - R1 / (D1 * n) * int(f2*sin(n*theta), theta, 0, 2*pi));
        
        An = 1 / (D1 * n * pi * R1^(n-1)) * (v1 + int(f2*cos(n*theta), theta, 0, 2*pi));
        Bn = 1 / (D1 * n * pi * R1^(n-1)) * (v2 + int(f2*sin(n*theta), theta, 0, 2*pi));

        C0 = -R1*log(R2) / (2*pi*D2) * v0;
        D0 = R1 / (2*pi*D2) * v0;

        Cn = 1 / (D2*n*pi*(R1^(n-1) + R2^(2*n)*R1^(-n-1))) * v1;
        Dn = -R2^(2*n) / (D2*n*pi*(R1^(n-1) + R2^(2*n)*R1^(-n-1))) * v1;

        Fn = 1 / (D2*n*pi*(R1^(n-1) + R2^(2*n)*R1^(-n-1))) * v2;
        Gn = -R2^(2*n) / (D2*n*pi*(R1^(n-1) + R2^(2*n)*R1^(-n-1))) * v2;
        
        STfunc{order + 1, 1} = 0;
        STfunc{order + 1, 2} = C0 + D0*log(r);
        
        for i = 1:maxfourier
            STfunc{order + 1, 1} = STfunc{order + 1, 1} + limit(An, n, i)*r^i*cos(i*theta) + limit(Bn, n, i)*r^i*sin(i*theta);
            STfunc{order + 1, 2} = STfunc{order + 1, 2} + (limit(Cn, n, i)*r^i + limit(Dn, n, i)*r^(-i))*cos(i*theta) + (limit(Fn, n, i)*r^i + limit(Gn, n, i)*r^(-i))*sin(i*theta);
        end
        
        STfunc{order + 1, 1} = simplify(STfunc{order + 1, 1}, 3);
        STfunc{order + 1, 2} = simplify(STfunc{order + 1, 2}, 3);
    end