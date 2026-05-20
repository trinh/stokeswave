clear
close all

N = 100; % Truncation
M = 30; % Epsilon iterations
kap = 0;
grav_coeff = 1;

eps = linspace(0.01,0.5, M);

p=0.25;
q = 1;
path = [p, p+q*1i, -p+q*1i, -p, p, p+q*1i, -p+q*1i, -p, p];

for k=1:length(eps)
    
% Initial guess - last two are mu and B (in that order!)
x0 = zeros(1,N-1);
x0(N) = 1;
x0(N+1) = 0.1;
ep = eps(k);

    
% Solve the ODE
options = optimset('Display', 'iter');
fwd = @(x)series_bern(x, ep, N, kap, grav_coeff);
[est,fval] = fsolve(fwd,x0,options);

% Save all the data
solseries = [];
solseries.N = N;
solseries.an = est(1:N-1);
solseries.mu = est(N);
solseries.B = est(N+1);
solseries.ep = ep;
solseries.kap = kap;
solseries.ycrest = sum(est(1:N-1)./(2*pi*(1:N-1)));


[f,g,gp,~,singloc] = gbern_contour(path,solseries);
location(k) = singloc;

[~,~,~,~,coeff(k)] = gbern_contour_coefficient(path,solseries,location(k));
%coeff
%abs(coeff)

end

hold on
plot(eps, imag(location))

figure(2);
plot(eps, abs(coeff));
