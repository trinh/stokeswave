% Get the coefficients of z'(f) for steep waves
% -----------------------------
% Written July 26 2015
% From Vanden-Broeck, 1986 Phys. Fluids

%clear
close all

savesol = 'false';

N = 120;
ep = sqrt(0.99); %ep^2 = 0.99

F = load('../../stokeswavedata/series_steep_N100_ep0.95.mat');
b = F.steepsolseries.bn;
mu = F.steepsolseries.mu;


usebeta = 1; % 1 to use beta formulation

% initial guess
x0 = zeros(1,N);
if usebeta == 0
    x0(1:N) = b(1:N);
    x0(N+1) = mu;
elseif usebeta == 1
    %x0(1:N-1) = b(1:N-1);
    x0 = zeros(1,N);
    x0(1) = 0.04119;
    x0(2) = 0.01252;
    x0(3) = 0.00606;
    x0(4) = 0.00360;
    x0(5) = 0.00240;
    x0(6) = 0.00173;
    x0(7) = 0.00131;
    x0(8) = 0.00102;
    x0(9) = 0.00083;
    x0(10) = 0.0068;
    x0(11) = 0.0057;
    x0(N) = 0.99; % This is beta
    x0(N+1) = mu;
end

%x0=est;

mytol = 1e-13;
options = optimset('Display', 'iter', 'TolFun', mytol, 'TolX', mytol);
fwd = @(x)vb_series_bern_steep(x, ep, N, usebeta);
[est,fval] = fsolve(fwd,x0,options); % Call solver

steepsolseries = [];
steepsolseries.N = N;
if usebeta == 1
    steepsolseries.bn = est(1:N-1);
    steepsolseries.beta = est(N);
elseif usebeta == 0
    steepsolseries.bn = est(1:N);
    steepsolseries.beta = 0;
end
steepsolseries.mu = est(N+1);
steepsolseries.ep = ep;

steepsolseries

%% Plot the solution

[phi, z, zp] = get_profile_steepsol(steepsolseries);

figure(1);
plot(real(z), imag(z), 'b.-'); hold all

%% Save the solution

switch savesol
    case 'true'
        name = ['../../stokeswavedata/series_steep_N', num2str(N), '_ep', num2str(ep), '.mat'];
        save(name, 'steepsolseries');
end


