function [f, x, y, xp, yp] = bdint_pathsolve_loadmesh(s, solbdint, solseries, meshdatbd)
% Solve ODE along a path 
% s(1) ---- s(2) ---- ... --- s(n-1) --- s(n)

f = [];
x = [];
y = [];
xp = [];
yp = [];

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
    fwd = @(s,U)bdint_bdfunc_loadmesh(gam(s), U, solbdint, solseries, meshdatbd)*dgam;
    
    if j == 1
        Z0 = [interp1(solbdint.phi, solbdint.x, real(z0)), ...
             interp1(solbdint.phi, solbdint.y, real(z0))];
    else 
        Z0 = [x(end), y(end)];
    end
    [ss, uu] = ode45(fwd, [0, 1], Z0);

    Fp = bdint_bdfunc_loadmesh(gam(ss).', uu.', solbdint, solseries, meshdatbd);    
    f = [f; gam(ss)];
    x = [x; uu(:,1)];
    y = [y; uu(:,2)];   
    xp = [xp; Fp(1,:).'];
    yp = [yp; Fp(2,:).'];
end
