function I = pvi_nyiri(f, fp, h, hp, hpp, x0, a, b)
% Principal value of int_{x_0-a}^{x_0+a} f(x)/[h(x) - h(x0)] dx

% Function to integrate. 
% We just "symmetrize" the integrand, to get rid of the singularity,
% but, to be able to integrate it numerically, 
% we have to provide a value at the singularity.
% Details: http://mat76.mat.uni-miskolc.hu/~mnotes/downloader.php?article=mmn_16.pdf

myint = @(f, a1, a2)integral(f, a1, a2);
    % I = quad(@gfunc, 0, a, 1e-12);
    % I = quadgk(@gfunc, 0, a);
    % I = rombquad(@gfunc, 0, a, 12);
    % /I = gaussquad(@gfunc, 0, L);

if (isreal(x0) == 1) && (x0 > a) && (x0 < b)
    % Estimate the derivatives of f and h
    % epsilon = 1e-12;
    % f1 = (f(x0 + epsilon) - f(x0-epsilon))/(2*epsilon);
    % h1 = (h(x0+epsilon) - h(x0-epsilon))/(2*epsilon);
    % h2 = (h(x0+epsilon) + h(x0-epsilon) - 2*h(x0))/epsilon^2;
    f1 = fp(x0);
    h1 = hp(x0);
    h2 = hpp(x0);    

    % Integrate in three parts
    % Take default value of L
    L = min(abs(b - x0), abs(a - x0));
    L = L/4;

    % First part from a to x0 - L
    I1 = myint(@(t)f(t)./(h(t) - h(x0)), a, x0 - L);

    % Second part from x0 - L to x0 + L
    I2 = myint(@gfunc, 0, L);

    % Third part from x0 + L to b
    I3 = myint(@(t)f(t)./(h(t) - h(x0)), x0 + L, b);

    I = I1 + I2 + I3;
else
    I = myint(@(t) f(t)./h(t), a, b);
end

% keyboard
    function g = gfunc(u)
        g = 0*u;
        g(u == 0) = (2 * f1 / h1 - f(x0) * h2 / h1^2);
        
        ind = u ~= 0;
        g(ind) = (f( x0 - u(ind) ) ./ ( h(x0 - u(ind)) - h(x0) ) + f( x0 + u(ind) ) ./ ( h(x0+u(ind)) - h(x0) ));
    end
end