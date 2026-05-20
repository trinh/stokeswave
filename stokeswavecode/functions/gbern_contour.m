function [f, g, gp, zp, singloc] = gbern_contour(s, solseries, g0)
% Solve ODE along a path
% s(1) ---- s(2) ---- ... --- s(n-1) --- s(n)

mytol = 1e-11;
options = odeset('RelTol', mytol, 'AbsTol', mytol);

if nargin < 3
    g0 = zfree(s(1), solseries) + zfree(-s(1), solseries);    
end

f = [];
g = [];
p_UHP = [];
gp = [];
integralval = 0;
for j = 1:length(s)-1
    
    z0 = s(j);
    z1 = s(j+1);
    
    % Compute straight line
    gam = @(s) z0 + s*(z1 - z0);
    dgam = (z1 - z0);
    
    % Now solve the ODE with initial conditions
    fwd = @(s,g)gbern_lower2upper(gam(s), g, solseries)*dgam;
    if j > 1
        g0 = g(end);
    end
    
%     [ss, gg] = ode45(fwd, [0, 1], g0, options);        
    [ss, gg] = ode113(fwd, [0, 1], g0, options);        
    
    % Also useful to get gp and p values
    [ggp, pp] = gbern_lower2upper(gam(ss), gg, solseries);
    
    p_UHP = [p_UHP; pp];
    gp = [gp; ggp];
    f = [f; gam(ss)];
    g = [g; gg];
    
   
    gg2 = gg - 2i*solseries.B*solseries.mu/(2*pi);
    integrand = dgam*gam(ss).*ggp./gg2;
    integralval = integralval + trapz(ss,integrand);

end


singloc = (1/(2*pi*1i))*(integralval);

zp = gp + p_UHP;
