function g = gbern_upper2lower_path
% GBERN_UPPER2LOWER_PATH

% =====================================
%         Load in the solution
F = load('data/series_N200_ep0.05.mat');
solseries = F.solseries;

F = load('data/meshdat_M100_ep0.05.mat');
meshdat = F.meshdat;
% =====================================

phi0 = -0.5;
s = [phi0, phi0 - 1i];

% Initial condition on g
phi = meshdat.UHP(1).f(1,:);
gaxis1 = meshdat.UHP(1).g(1,:);
gaxis2 = meshdat.UHP(2).g(1,:);
g0 = interp1(phi, gaxis2, phi0);

%%
f = [];
g = [];
for j = 1:length(s)-1
    
    z0 = s(j);
    z1 = s(j+1);
    
    % Compute straight line
    gam = @(s) z0 + s*(z1 - z0);
    dgam = (z1 - z0);
    
    % Now solve the ODE with initial conditions
    fwd = @(s,g)gbern_upper2lower_pathode(gam(s), g, solseries, phi, gaxis1)*dgam;
    if j > 1
        g0 = g(end);
    end
   
    [ss, gg] = ode45(fwd, [0, 1], g0);
    f = [f; gam(ss)];
    g = [g; gg];
    
    keyboard
end
