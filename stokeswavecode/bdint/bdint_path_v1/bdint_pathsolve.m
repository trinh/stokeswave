function [f, x, y] = bdint_pathsolve(s, solbdint, solseries)
% Solve ODE along a path 
% s(1) ---- s(2) ---- ... --- s(n-1) --- s(n)

f = [];
x = [];
y = [];

mytol = 1e-6;
options = odeset('RelTol', mytol, 'AbsTol', mytol);
for j = 1:length(s)-1
    
    z0 = s(j);
    z1 = s(j+1);
                  
    if (isnan(z0) == 1) || (isnan(z1) == 1)        
        continue;
    end
    
    % Compute straight line
    gam = @(s) z0 + s*(z1 - z0);
    dgam = (z1 - z0);
    
    % Now solve the ODE with initial conditions
    fwd = @(s,U)bdfunc_complex(gam(s), U, solbdint, solseries)*dgam;
    
    if j == 1
        Z0 = [interp1(solbdint.phi, solbdint.x, real(z0)), ...
             interp1(solbdint.phi, solbdint.y, real(z0))];
    else 
        Z0 = [x(end), y(end)];
    end
    [ss, uu] = ode45(fwd, [0, 1], Z0, options);    
    
    f = [f; gam(ss)];
    x = [x; uu(:,1)];
    y = [y; uu(:,2)];        
end
