% Get the coefficients of z'(f) via series truncation and collocation
% -----------------------------
% Written June-July 2015

clear
close all

% =====================================
% Load in the solution for an initial guess
%F = load('data/series_N100_ep0.5.mat');
%solseriesguess = F.solseries;
% =====================================


N = 100; % N-1 is no. of terms in power series
% ep = 0.7746; % JMVB 1986 (check mu ~ 1.22)
ep = 0.1;
kap = 0;
grav_coeff = 1;

% Initial guess - last two are mu and B (in that order!)
x0 = zeros(1,N-1);
x0(N) = 1;
x0(N+1) = 0.1;


% Or take from loaded solution
% x0(1:N-1) = solseriesguess.an;
% x0(N) = solseriesguess.mu;
% x0(N+1) = solseriesguess.B;



% Solve the ODE
mytol = 1e-13;
options = optimset('Display', 'iter', 'TolFun', mytol, 'TolX', mytol);
fwd = @(x)series_bern(x, ep, N, kap, grav_coeff);
[est,fval] = fsolve(fwd,x0,options); % Call solver


%%
% Save all the data
solseries = [];
solseries.N = N;
solseries.an = est(1:N-1);
solseries.mu = est(N);
solseries.B = est(N+1);
solseries.ep = ep;
solseries.kap = kap;
solseries.ycrest = sum(est(1:N-1)./(2*pi*(1:N-1)));
solseries.H = imag(zfree(0,solseries)-zfree(0.5,solseries));

name = ['../../stokeswavedata/series_N', num2str(N), '_ep', num2str(ep), '.mat'];
save(name, 'solseries');


%% Plot the solution

phi = linspace(-0.5,0.5,100);
figure(1);
plot(real(zfree(phi,solseries)), imag(zfree(phi,solseries)));
axis equal


