function [t, path] = het_disc2(R1star, R2, P, v_0, delta, tau, n, do_vis)

    P_vec = [P, P(end)];
    t = 0;
    v = v_0;
    if(do_vis)
        path = [v(1), v(2), 1];
    else
        path = [];
    end

    [theta, r] = cart2pol(v(1), v(2));
    
    while (r < R2)
        
        [new_x, new_y, p, new_theta, new_r] = get_new_pos(R1star, R2, P_vec, v(1), v(2), delta, n);
        prob = rand;
        
        %i = find((p - prob) > 0, 1) - 1;
        
        a = 0;
        b = length(p);
        while(a < b)
            i = floor((a + b) / 2);
            if(p(i + 1) < prob)
                a = i + 1;
            else
                b = i;
            end
        end
        
        if(a == length(p) - 1)
            t = t + tau;
            if(do_vis)
                path(end, 3) = path(end, 3) + 1;
            end
        else
            v = [new_x(a), new_y(a)];
            theta = new_theta(a);
            r = new_r(a);
            t = t + tau;
            if(do_vis)
                path = [path; v(1), v(2), 1];
            end
        end
                
    end

end

function [x2, y2, p, theta2, r2] = get_new_pos(R1star,R2, P_vec, x1, y1, delta, n)
    
    theta1 = 2 * pi * (0:n - 1) / n;

    x2 = x1 + delta * cos(theta1);
    y2 = y1 + delta * sin(theta1);
    
    midx = x1 + 0.5 * delta * cos(theta1);
    midy = y1 + 0.5 * delta * sin(theta1);
    
    [theta2, r2] = cart2pol(x2, y2);
    [midtheta, midr] = cart2pol(midx, midy);
    
    p = zeros(1, length(theta1) + 1);
    for i = 1:length(p) - 1
%         for j = 1:length(R_vec) - 1
%             if R_vec{j}(midtheta(i)) < midr(i) && midr(i) < R_vec{j+1}(midtheta(i))
%                 p(i + 1) = p(i) + P_vec(j);
%                 break
%             end
            if midr(i) < R1star(midtheta(i))
                p(i + 1) = p(i) + P_vec(1);
            else
                p(i + 1) = p(i) + P_vec(2);
            end
%         end
    end    
    p = [p / n, 1];
end
