
p=0.25;
q = 1;
path = [p, p+q*1i, -p+q*1i, -p, p, p+q*1i, -p+q*1i, -p, p];

ep = 0.3158;
kap = 0;
grav_coeff = 1;

% Initial guess - last two are mu and B (in that order!)

x0 = zeros(1,N-1);
x0(N) = 1;
x0(N+1) = 0.1;

% Solve the ODE
options = optimset('Display', 'iter');
fwd = @(x)series_bern(x, ep, N, kap, grav_coeff);
[est,fval] = fsolve(fwd,x0,options); % Call solver

% Get the wave height
solseries = [];
solseries.N = N;
solseries.an = est(1:N-1);
solseries.mu = est(N);
solseries.B = est(N+1);
solseries.ep = ep;
solseries.kap = kap;
solseries.ycrest = sum(est(1:N-1)./(2*pi*(1:N-1)));

[~,~,~,~,singloc] = gbern_contour(path,solseries);
[~,~,~,~,coeff] = gbern_contour_coefficient(path,solseries,singloc);