function [f, G, integral] = odecontour(s,singloc)
% Solve ODE along a path 
% s(1) ---- s(2) ---- ... --- s(n-1) --- s(n)

    f = [];
    G = [];
    integral = 0;

    for j = 1:length(s)-1
        f0 = s(j);
        f1 = s(j+1);

        % Compute straight line
        gam = @(t) f0 + t*(f1 - f0);
        dgam = (f1 - f0);

        % Now calculate dw/dz
        fwd = @(t,y)myode(gam(t), y)*dgam;
        if j == 1
            g0 = (s(1)-singloc)^(1/2);
        else 
            g0 = G(end);
        end

        [tt, gg] = ode45(fwd, [0, 1], g0);

        ggp = myode(gam(tt), gg);

        f = [f; gam(tt)];
        G = [G; gg];

        integrand = dgam*gam(tt).*ggp./gg;
        integral = integral + trapz(tt,integrand);
    end

end

function gp = myode(f, G)
   gp = 0.5./G;
end

    

