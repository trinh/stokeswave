% Relationship between wave height and mu vs. ep
clear all 
close all

N = 25; % Truncation
M = 50; %ep iterations

eps = linspace(0.6,0.9,M);
est = zeros(N+1,M);

for k=1:M
ep = eps(k)
kap = 0;
grav_coeff = 1;

% Initial guess - last two are mu and B (in that order!)
if k == 1
    x0 = zeros(1,N-1);
    x0(N) = 1;
    x0(N+1) = 0.1;
else
    x0 = est(:,k-1);
    x0 = transpose(x0);
end

% Solve the ODE
mytol = 1e-13;
options = optimset('Display', 'iter', 'TolFun', mytol, 'TolX', mytol);
fwd = @(x)series_bern(x, ep, N, kap, grav_coeff);
[est(:,k),fval] = fsolve(fwd,x0,options); % Call solver

% Get the wave height
solseries = [];
solseries.N = N;
solseries.an = est(1:N-1,k);
solseries.mu = est(N,k);
solseries.B = est(N+1,k);
solseries.ep = ep;
solseries.kap = kap;
solseries.ycrest = sum(est(1:N-1)./(2*pi*(1:N-1)));
height(k) = imag(zfree(0,solseries)-zfree(0.5,solseries));

end

%%
% Plot ep against mu
 
figure(1);
plot(eps,est(N,:));
title('\mu')

figure(2);
plot(eps,height);
title('Height')
