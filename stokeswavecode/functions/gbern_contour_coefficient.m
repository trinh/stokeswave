function [f, g, gp, zp, coeff] = gbern_contour_coefficient(s, solseries, singloc, g0)
% Solve ODE along a path
% s(1) ---- s(2) ---- ... --- s(n-1) --- s(n)
% Much the same as gbern_contour except here we calculate the coefficient
% integral instead. Need the singularity location as an input

if nargin < 4
    g0 = zfree(s(1), solseries) + zfree(-s(1), solseries);    
end

f = [];
g = [];
p_UHP = [];
gp = [];
integralval = 0;
c_intval = 0;
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
    
    [ss, gg] = ode45(fwd, [0, 1], g0);        
    
%     gg = gg - 2*solseries.ycrest;
    % Also useful to get gp and p values
    [ggp, pp] = gbern_lower2upper(gam(ss), gg, solseries);
    
    p_UHP = [p_UHP; pp];
    gp = [gp; ggp];
    f = [f; gam(ss)];
    g = [g; gg];
    
    zz = gg - zfree(-gam(ss), solseries);
    zzp = ggp + pp;
   
    gg2 = gg - 2i*solseries.B*solseries.mu/(2*pi);
    integrand = dgam*gam(ss).*ggp./gg2;
    integralval = integralval + trapz(ss,integrand);
    
    c_integrand = ((gam(ss)-singloc).^(-3/2)).*gg2*dgam;
    c_intval = c_intval + trapz(ss,c_integrand);

end


coeff = (1/(4*pi*1i))*(c_intval);

zp = gp + p_UHP;
