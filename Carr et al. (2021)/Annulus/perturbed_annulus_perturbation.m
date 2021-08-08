function T = perturbed_annulus_perturbation(R_inner, R_outer, g_inner, g_outer, ...
    D, perturb_degree, N, epsilon, inner_type, backup_inner, backup_outer, backup_D)
%AnnulusPerturbation: Finds the perturbation solution T = T0 + epsilon T1 +
%epsilon^2 T2 + ... + epsilon^n Tn to the problem of mean first passage
%time in a perturbed annulus with outer boundary described by
%R_outer^*(theta) = R_outer(1 + epsilon g_outer(theta)) and the inner
%boundary by R_inner^*(theta) = R_inner(1 + epsilon g_inner(theta)). Here,
%R_outer is the unperturbed outer radius of the innulus and R_inner is the
%unperturbed inner radius. g_inner and g_outer ad smooth 2pi periodic
%functions. The diffusivity is D, we take N terms in the Fourier series,
%the perturbation parameter is epsilon << 1, and the number of terms in the
%perturbation series is n. Inner_type is 'reflect' if reflective, and
%'absorb', if we absorb at the inner boundary.
syms r; assume(r > 0);
syms theta;
syms n; assume(n, 'integer'); assumeAlso(n > 0);
m = perturb_degree + 1;
syms T [1 m]
syms Hk_outer(theta)
syms Hk_inner(theta)
temp_err = warning('error', 'symbolic:solve:warnmsg3');
if strcmpi(inner_type, 'reflect')
    T(1) = (R_outer^2 - r.^2 + 2*R_inner^2 * log(r/R_outer))/(4*D);
    syms dg(theta)
    dg(theta) = diff(g_inner, theta);
    if perturb_degree >= 1
        Hk_outer(theta) = - R_outer*g_outer(theta).*subs(diff(T(1), r), r, R_outer);
        Hk_inner(theta) = -R_inner*g_inner(theta).*subs(diff(T(1), r, 2), r, R_inner) - 2*g_inner(theta) .* subs(diff(T(1), r), r, R_inner);
        C01 = 1/(2*D*sym(pi))*intfix((R_outer^2-R_inner^2)*g_outer(theta)-2*R_inner^2*log(R_outer)*g_inner(theta), theta, 0, 2*sym(pi));
        D01 = R_inner^2/(D*sym(pi))*intfix(g_inner(theta),theta,0,2*sym(pi));
        Vn = R_inner^(-n-1)*R_outer^n+R_inner^(n-1)*R_outer^(-n);
        An1 = 1/(2*D*sym(pi)*n*Vn)*intfix((n*R_inner^(-n-1)*(R_outer^2-R_inner^2)*g_outer(theta)+2*R_outer^(-n)*R_inner*g_inner(theta))*sin(n*theta),theta,0,2*sym(pi));
        Bn1 = 1/(2*D*sym(pi)*n*Vn)*intfix((n*R_inner^(n-1)*(R_outer^2-R_inner^2)*g_outer(theta)-2*R_outer^(n)*R_inner*g_inner(theta))*sin(n*theta),theta,0,2*sym(pi));
        Cn1 = 1/(2*D*sym(pi)*n*Vn)*intfix((n*R_inner^(-n-1)*(R_outer^2-R_inner^2)*g_outer(theta)+2*R_outer^(-n)*R_inner*g_inner(theta))*cos(n*theta),theta,0,2*sym(pi));
        Dn1 = 1/(2*D*sym(pi)*n*Vn)*intfix((n*R_inner^(n-1)*(R_outer^2-R_inner^2)*g_outer(theta)-2*R_outer^(n)*R_inner*g_inner(theta))*cos(n*theta),theta,0,2*sym(pi));
        [~, dAn1] = numden(An1);
        [~, dBn1] = numden(Bn1);
        [~, dCn1] = numden(Cn1);
        [~, dDn1] = numden(Dn1);
        try
            d1 = solve(dAn1 == 0, n);
        catch
            warning('Failed to solve for roots, setting inner_radius and outer_radius to %g and %g and trying again.', backup_inner, backup_outer)
            try
                d1 = solve(subs(dAn1, [R_inner R_outer], [backup_inner backup_outer]) == 0, n);
            catch
                warning('Second attempt failed, also setting D to %g and trying again.', backup_D)
                d1 = solve(subs(dAn1, [R_inner R_outer D], [backup_inner backup_outer backup_D]) == 0, n);
            end
        end
        try
            d2 = solve(dBn1 == 0, n);
        catch
            warning('Failed to solve for roots, setting inner_radius and outer_radius to %g and %g and trying again.', backup_inner, backup_outer)
            try
                d2 = solve(subs(dBn1, [R_inner R_outer], [backup_inner backup_outer]) == 0, n);
            catch
                warning('Second attempt failed, also setting D to %g and trying again.', backup_D)
                d2 = solve(subs(dBn1, [R_inner R_outer D], [backup_inner backup_outer backup_D]) == 0, n);
            end
        end
        try
            d3 = solve(dCn1 == 0, n);
        catch
            warning('Failed to solve for roots, setting inner_radius and outer_radius to %g and %g and trying again.', backup_inner, backup_outer)
            try
                d3 = solve(subs(dCn1, [R_inner R_outer], [backup_inner backup_outer]) == 0, n);
            catch
                warning('Second attempt failed, also setting D to %g and trying again.', backup_D)
                d3 = solve(subs(dCn1, [R_inner R_outer D], [backup_inner backup_outer backup_D]) == 0, n);
            end
        end
        try
            d4 = solve(dDn1 == 0, n);
        catch
            warning('Failed to solve for roots, setting inner_radius and outer_radius to %g and %g and trying again.', backup_inner, backup_outer)
            try
                d4 = solve(subs(dDn1, [R_inner R_outer], [backup_inner backup_outer]) == 0, n);
            catch
                warning('Second attempt failed, also setting D to %g and trying again.', backup_D)
                d4 = solve(subs(dDn1, [R_inner R_outer D], [backup_inner backup_outer backup_D]) == 0, n);
            end
        end
        divideZeroIdx = double(union(union(round(d1), round(d2)), union(round(d3), round(d4))))';
        divideZeroIdx = divideZeroIdx(divideZeroIdx > 0 & divideZeroIdx <= N);
        T(2) = 1/2 * (C01 + D01 * log(r));
        if ~isempty(divideZeroIdx)
            for nn = divideZeroIdx
                Ann1 = 1/(2*D*sym(pi)*nn*subs(Vn, n, nn))*intfix((nn*R_inner^(-nn-1)*(R_outer^2-R_inner^2)*g_outer(theta)+2*R_outer^(-nn)*R_inner*g_inner(theta))*sin(nn*theta),theta,0,2*sym(pi));
                Bnn1 = 1/(2*D*sym(pi)*nn*subs(Vn, n, nn))*intfix((nn*R_inner^(nn-1)*(R_outer^2-R_inner^2)*g_outer(theta)-2*R_outer^(nn)*R_inner*g_inner(theta))*sin(nn*theta),theta,0,2*sym(pi));
                Cnn1 = 1/(2*D*sym(pi)*nn*subs(Vn, n, nn))*intfix((nn*R_inner^(-nn-1)*(R_outer^2-R_inner^2)*g_outer(theta)+2*R_outer^(-nn)*R_inner*g_inner(theta))*cos(nn*theta),theta,0,2*sym(pi));
                Dnn1 = 1/(2*D*sym(pi)*nn*subs(Vn, n, nn))*intfix((nn*R_inner^(nn-1)*(R_outer^2-R_inner^2)*g_outer(theta)-2*R_outer^(nn)*R_inner*g_inner(theta))*cos(nn*theta),theta,0,2*sym(pi));
                T(2) = T(2) + (Cnn1*r^nn + Dnn1*r^(-nn))*cos(nn*theta) + (Ann1*r^nn + Bnn1*r^(-nn))*sin(nn*theta);
            end
        end
        for nn = setdiff(1:N, divideZeroIdx)
            Cnn1 = subs(Cn1, n, nn);
            Dnn1 = subs(Dn1, n, nn);
            Bnn1 = subs(Bn1, n, nn);
            Ann1 = subs(An1, n, nn);
            T(2) = T(2) + (Cnn1*r^nn + Dnn1*r^(-nn))*cos(nn*theta) + (Ann1*r^nn + Bnn1*r^(-nn))*sin(nn*theta);
        end
        if perturb_degree >= 2
            for k = 2:perturb_degree
                Hk_outer(theta) = 0;
                for j=1:k
                    Hk_outer(theta) = Hk_outer(theta) - (R_outer*g_outer(theta))^j/factorial(j)*subs(diff(T(k-j+1), r, j), r, R_outer);
                end
                Hk_inner(theta) = -R_inner^k*g_inner(theta).^k/factorial(k)*subs(diff(T(1), r, k+1), r, R_inner);
                Hk_inner(theta) = Hk_inner(theta) - R_inner^(k-1)*g_inner(theta).^(k-1)/factorial(k-1)*subs(diff(T(2), r, k), r, R_inner);
                Hk_inner(theta) = Hk_inner(theta) - 2*R_inner^(k-1)*g_inner(theta).^k/factorial(k-1)*subs(diff(T(1), r, k), r, R_inner);
                for j=1:k-2
                    Hk_inner(theta) = Hk_inner(theta) - R_inner^j*g_inner(theta).^j/factorial(j)*subs(diff(T(k-j+1), r, j+1), r, R_inner);
                end
                for j=0:k-2
                    Hk_inner(theta) = Hk_inner(theta) + R_inner^(j-1)*g_inner(theta).^j.*dg(theta)/factorial(j)*subs(diff(diff(T(k-j), theta), r, j), r, R_inner);
                    Hk_inner(theta) = Hk_inner(theta) - 2*R_inner^j*g_inner(theta).^(j+1)/factorial(j)*subs(diff(T(k-j), r, j+1), r, R_inner);
                    Hk_inner(theta) = Hk_inner(theta) - R_inner^j*g_inner(theta).^(j+2)/factorial(j)*subs(diff(T(k-j-1), r, j+1), r, R_inner);
                end
                C0k = 1/sym(pi) * intfix(Hk_outer(theta) - R_inner*log(R_outer)*Hk_inner(theta), theta, 0, 2*sym(pi));
                D0k = 1/sym(pi) * intfix(R_inner*Hk_inner(theta), theta, 0, 2*sym(pi));
                Vn = R_inner^(-n-1)*R_outer^n + R_inner^(n-1)*R_outer^(-n);
                Ank = 1/(n*sym(pi)*Vn) * intfix((n*R_inner^(-n-1)*Hk_outer(theta) + R_outer^(-n)*Hk_inner(theta))*sin(n*theta), theta, 0, 2*sym(pi));
                Bnk = 1/(n*sym(pi)*Vn) * intfix((n*R_inner^(n-1)*Hk_outer(theta) - R_outer^(-n)*Hk_inner(theta))*sin(n*theta), theta, 0, 2*sym(pi));
                Cnk = 1/(n*sym(pi)*Vn) * intfix((n*R_inner^(-n-1)*Hk_outer(theta) + R_outer^(-n)*Hk_inner(theta))*cos(n*theta), theta, 0, 2*sym(pi));
                Dnk = 1/(n*sym(pi)*Vn) * intfix((n*R_inner^(n-1)*Hk_outer(theta) - R_outer^(-n)*Hk_inner(theta))*cos(n*theta), theta, 0, 2*sym(pi));
                [~, dAnk] = numden(Ank);
                [~, dBnk] = numden(Bnk);
                [~, dCnk] = numden(Cnk);
                [~, dDnk] = numden(Dnk);
                try
                    d1 = solve(dAnk == 0, n);
                catch
                    warning('Failed to solve for roots, setting inner_radius and outer_radius to %g and %g and trying again.', backup_inner, backup_outer)
                    try
                        d1 = solve(subs(dAnk, [R_inner R_outer], [backup_inner backup_outer]) == 0, n);
                    catch
                        warning('Second attempt failed, also setting D to %g and trying again.', backup_D)
                        d1 = solve(subs(dAnk, [R_inner R_outer D], [backup_inner backup_outer backup_D]) == 0, n);
                    end
                end
                try
                    d2 = solve(dBnk == 0, n);
                catch
                    warning('Failed to solve for roots, setting inner_radius and outer_radius to %g and %g and trying again.', backup_inner, backup_outer)
                    try
                        d2 = solve(subs(dBnk, [R_inner R_outer], [backup_inner backup_outer]) == 0, n);
                    catch
                        warning('Second attempt failed, also setting D to %g and trying again.', backup_D)
                        d2 = solve(subs(dBnk, [R_inner R_outer D], [backup_inner backup_outer backup_D]) == 0, n);
                    end
                end
                try
                    d3 = solve(dCnk == 0, n);
                catch
                    warning('Failed to solve for roots, setting inner_radius and outer_radius to %g and %g and trying again.', backup_inner, backup_outer)
                    try
                        d3 = solve(subs(dCnk, [R_inner R_outer], [backup_inner backup_outer]) == 0, n);
                    catch
                        warning('Second attempt failed, also setting D to %g and trying again.', backup_D)
                        d3 = solve(subs(dCnk, [R_inner R_outer D], [backup_inner backup_outer backup_D]) == 0, n);
                    end
                end
                try
                    d4 = solve(dDnk == 0, n);
                catch
                    warning('Failed to solve for roots, setting inner_radius and outer_radius to %g and %g and trying again.', backup_inner, backup_outer)
                    try
                        d4 = solve(subs(dDnk, [R_inner R_outer], [backup_inner backup_outer]) == 0, n);
                    catch
                        warning('Second attempt failed, also setting D to %g and trying again.', backup_D)
                        d4 = solve(subs(dDnk, [R_inner R_outer D], [backup_inner backup_outer backup_D]) == 0, n);
                    end
                end
                divideZeroIdx = double(union(union(round(d1), round(d2)), union(round(d3), round(d4))))';
                divideZeroIdx = divideZeroIdx(divideZeroIdx > 0 & divideZeroIdx <= N);
                T(k+1) = 1/2 * (C0k + D0k * log(r));
                if ~isempty(divideZeroIdx)
                    for nn = divideZeroIdx
                        Annk = 1/(nn*sym(pi)*subs(Vn, n,nn)) * intfix((nn*R_inner^(-nn-1)*Hk_outer(theta) + R_outer^(-nn)*Hk_inner(theta))*sin(nn*theta), theta, 0, 2*sym(pi));
                        Bnnk = 1/(nn*sym(pi)*subs(Vn,n, nn)) * intfix((nn*R_inner^(nn-1)*Hk_outer(theta) - R_outer^(-nn)*Hk_inner(theta))*sin(nn*theta), theta, 0, 2*sym(pi));
                        Cnnk = 1/(nn*sym(pi)*subs(Vn,n, nn)) * intfix((nn*R_inner^(-nn-1)*Hk_outer(theta) + R_outer^(-nn)*Hk_inner(theta))*cos(nn*theta), theta, 0, 2*sym(pi));
                        Dnnk = 1/(nn*sym(pi)*subs(Vn,n, nn)) * intfix((nn*R_inner^(nn-1)*Hk_outer(theta) - R_outer^(-nn)*Hk_inner(theta))*cos(nn*theta), theta, 0, 2*sym(pi));
                        T(k+1) = T(k+1) + (Cnnk*r^nn + Dnnk*r^(-nn))*cos(nn*theta) + (Annk*r^nn + Bnnk*r^(-nn))*sin(nn*theta);
                    end
                end
                for nn = setdiff(1:N, divideZeroIdx)
                    Cnnk = subs(Cnk, n, nn);
                    Dnnk = subs(Dnk, n, nn);
                    Bnnk = subs(Bnk, n, nn);
                    Annk = subs(Ank, n, nn);
                    T(k+1) = T(k+1) + (Cnnk*r^nn + Dnnk*r^(-nn))*cos(nn*theta) + (Annk*r^nn + Bnnk*r^(-nn))*sin(nn*theta);
                end
            end
        end
        T = sum(epsilon.^(0:perturb_degree) .* T);
    end
elseif strcmpi(inner_type, 'absorb')
    T(1) = 1/(4*D*log(R_inner/R_outer)) * (R_inner^2*log(r/R_outer) + R_outer^2*log(R_inner./r) + r^2*log(R_outer/R_inner));
    for k = 1:perturb_degree
        Hk_outer(theta) = 0;
        Hk_inner(theta) = 0;
        for j = 1:k
            Hk_outer(theta) = Hk_outer(theta) -(R_outer.*g_outer(theta)).^j/factorial(j)*subs(diff(T(k-j+1),r,j),r,R_outer);
            Hk_inner(theta) = Hk_inner(theta) - (R_inner*g_inner(theta)).^j/factorial(j)*subs(diff(T(k-j+1),r,j),r,R_inner);
        end
        C0 = 1/(sym(pi)*log(R_outer/R_inner))*intfix(Hk_inner(theta)*log(R_outer)-Hk_outer(theta)*log(R_inner),theta,0,2*sym(pi));
        D0 = 1/(sym(pi)*log(R_outer/R_inner))*intfix(Hk_outer(theta)-Hk_inner(theta), theta, 0, 2*sym(pi));
        An = 1/(sym(pi)*(R_inner^(2*n)-R_outer^(2*n)))*intfix((R_inner^n*Hk_inner(theta)-R_outer^n*Hk_outer(theta)).*sin(n*theta),theta,0,2*sym(pi));
        Bn = 1/(sym(pi)*(R_inner^(-2*n)-R_outer^(-2*n)))*intfix((R_inner^(-n)*Hk_inner(theta)-R_outer^(-n)*Hk_outer(theta)).*sin(n*theta),theta,0,2*sym(pi));
        Cn = 1/(sym(pi)*(R_inner^(2*n)-R_outer^(2*n)))*intfix((R_inner^n*Hk_inner(theta)-R_outer^n*Hk_outer(theta)).*cos(n*theta),theta,0,2*sym(pi));
        Dn = 1/(sym(pi)*(R_inner^(-2*n)-R_outer^(-2*n)))*intfix((R_inner^(-n)*Hk_inner(theta)-R_outer^(-n)*Hk_outer(theta)).*cos(n*theta),theta,0,2*sym(pi));
        [~, dAn] = numden(An);
        [~, dBn] = numden(Bn);
        [~, dCn] = numden(Cn);
        [~, dDn] = numden(Dn);
        try
            d1 = solve(dAn == 0, n);
        catch
            warning('Failed to solve for roots, setting inner_radius and outer_radius to %g and %g and trying again.', backup_inner, backup_outer)
            try
                d1 = solve(subs(dAn, [R_inner R_outer], [backup_inner backup_outer]) == 0, n);
            catch
                warning('Second attempt failed, also setting D to %g and trying again.', backup_D)
                d1 = solve(subs(dAn, [R_inner R_outer D], [backup_inner backup_outer backup_D]) == 0, n);
            end
        end
        try
            d2 = solve(dBn == 0, n);
        catch
            warning('Failed to solve for roots, setting inner_radius and outer_radius to %g and %g and trying again.', backup_inner, backup_outer)
            try
                d2 = solve(subs(dBn, [R_inner R_outer], [backup_inner backup_outer]) == 0, n);
            catch
                warning('Second attempt failed, also setting D to %g and trying again.', backup_D)
                d2 = solve(subs(dBn, [R_inner R_outer D], [backup_inner backup_outer backup_D]) == 0, n);
            end
        end
        try
            d3 = solve(dCn == 0, n);
        catch
            warning('Failed to solve for roots, setting inner_radius and outer_radius to %g and %g and trying again.', backup_inner, backup_outer)
            try
                d3 = solve(subs(dCn, [R_inner R_outer], [backup_inner backup_outer]) == 0, n);
            catch
                warning('Second attempt failed, also setting D to %g and trying again.', backup_D)
                d3 = solve(subs(dCn, [R_inner R_outer D], [backup_inner backup_outer backup_D]) == 0, n);
            end
        end
        try
            d4 = solve(dDn == 0, n);
        catch
            warning('Failed to solve for roots, setting inner_radius and outer_radius to %g and %g and trying again.', backup_inner, backup_outer)
            try
                d4 = solve(subs(dDn, [R_inner R_outer], [backup_inner backup_outer]) == 0, n);
            catch
                warning('Second attempt failed, also setting D to %g and trying again.', backup_D);
                d4 = solve(subs(dDn, [R_inner R_outer D], [backup_inner backup_outer backup_D]) == 0, n);
            end
        end
        divideZeroIdx = double(union(union(round(d1), round(d2)), union(round(d3), round(d4))))';
        divideZeroIdx = divideZeroIdx(divideZeroIdx > 0 & divideZeroIdx <= N);
        T(k+1) = 1/2 * (C0 + D0 * log(r));
        if ~isempty(divideZeroIdx)
            for nn = divideZeroIdx
                Ann = 1/(sym(pi)*(R_inner^(2*nn)-R_outer^(2*nn)))*intfix((R_inner^nn*Hk_inner(theta)-R_outer^nn*Hk_outer(theta)).*sin(nn*theta),theta,0,2*sym(pi));
                Bnn = 1/(sym(pi)*(R_inner^(-2*nn)-R_outer^(-2*nn)))*intfix((R_inner^(-nn)*Hk_inner(theta)-R_outer^(-nn)*Hk_outer(theta)).*sin(nn*theta),theta,0,2*sym(pi));
                Cnn = 1/(sym(pi)*(R_inner^(2*nn)-R_outer^(2*nn)))*intfix((R_inner^nn*Hk_inner(theta)-R_outer^nn*Hk_outer(theta)).*cos(nn*theta),theta,0,2*sym(pi));
                Dnn = 1/(sym(pi)*(R_inner^(-2*nn)-R_outer^(-2*nn)))*intfix((R_inner^(-nn)*Hk_inner(theta)-R_outer^(-nn)*Hk_outer(theta)).*cos(nn*theta),theta,0,2*sym(pi));
                T(k+1) = T(k+1) + (Cnn*r^nn + Dnn*r^(-nn))*cos(nn*theta) + (Ann*r^nn + Bnn*r^(-nn))*sin(nn*theta);
            end
        end
        for nn = setdiff(1:N, divideZeroIdx)
            Cnn = subs(Cn, n, nn);
            Dnn = subs(Dn, n, nn);
            Bnn = subs(Bn, n, nn);
            Ann = subs(An, n, nn);
            T(k+1) = T(k+1) + (Cnn*r^nn + Dnn*r^(-nn))*cos(nn*theta) + (Ann*r^nn + Bnn*r^(-nn))*sin(nn*theta);
        end
    end
    T = sum(epsilon.^(0:perturb_degree) .* T);
end
warning(temp_err);
function I = intfix(expr, var, a, b)
I = int(expr, var);
I = simplifyFraction(subs(I, var, b) - subs(I, var, a));
%end