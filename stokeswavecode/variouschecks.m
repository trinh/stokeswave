% Relationship between wave height and mu vs. ep
% Also now looks at the coefficient and location of the "nearest" singularity
clear
close all

N = 20; % Truncation
M = 20; % ep iterations

eps = linspace(0.1,0.9,M);
kap = 0;
grav_coeff = 1;
est = zeros(N+1,M);

p=0.25;
q = 1;
path = [p, p+q*1i, -p+q*1i, -p, p, p+q*1i, -p+q*1i, -p, p];

for k=1:M
ep = eps(k);

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
options = optimset('Display', 'iter');
fwd = @(x)series_bern(x, ep, N, kap, grav_coeff);
[est(:,k),fval] = fsolve(fwd,x0,options); % Call solver

% Get the wave height
solseries = [];
solseries.N = N;
solseries.an = transpose(est(1:N-1,k));
solseries.mu = est(N,k);
solseries.B = est(N+1,k);
solseries.ep = ep;
solseries.kap = kap;
solseries.ycrest = sum(transpose(est(1:N-1,k))./(2*pi*(1:N-1)));
height(k) = imag(zfree(0,solseries)-zfree(0.5,solseries));

[~,~,~,~,singloc(k)] = gbern_contour(path,solseries);
[~,~,~,~,coeff(k)] = gbern_contour_coefficient(path,solseries,singloc(k));

spectral_estimate(k) = (1/(2*pi))*((1/((N-1)))*((-0.5)*log((N-1)) - log(abs(est(N-1,k)))));

end


% Now we find epsilon vs. mu for the steep solver

F = load('../stokeswavedata/series_steep_N100_ep0.95.mat');
b = F.steepsolseries.bn;
mu = F.steepsolseries.mu;

eps2 = linspace(0.995, 0.7, M);
est2 = zeros(N+1,M);

for k2=1:M
usebeta = 1; % 1 to use beta formulation
% initial guess
x0_2 = zeros(1,N);    
x0_2(1) = 0.04119;
x0_2(2) = 0.01252;
x0_2(3) = 0.00606;
x0_2(4) = 0.00360;
x0_2(5) = 0.00240;
x0_2(6) = 0.00173;
x0_2(7) = 0.00131;
x0_2(8) = 0.00102;
x0_2(9) = 0.00083;
x0_2(10) = 0.0068;
x0_2(11) = 0.0057;
x0_2(N) = 0.99; % This is beta
x0_2(N+1) = 1.1;  % This is mu

if k2==1
else
    x0_2 = est2(:,k2-1);
end

options = optimset('Display', 'iter');
fwd = @(x)vb_series_bern_steep(x, eps2(k2), N, usebeta);
[est2(:,k2),fval] = fsolve(fwd,x0_2,options); % Call solver
end

%%
% Plot stuff
close all

figure(1);
plot(eps,est(N,:));
hold on
plot(eps2,est2(N+1,:));
title('\mu')
legend('Standard series solver', 'steep series');

figure(2);
plot(eps,height);
title('Height')

figure(3);
plot(eps,imag(singloc));
hold on
plot(eps,spectral_estimate);
title('Singularity Location');
legend('Contour integral', 'Spectral method');

figure(4);
plot(eps,abs(coeff));
hold on
plot(eps,real(coeff));
plot(eps,imag(coeff));
title('Singularity Coefficient');
legend('abs', 'real part', 'imaginary part');
