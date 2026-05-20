% Get the coefficients of z'(f) for steep waves - iterating through epsilon
% -----------------------------
%
%   Re-written by Phil Sept 2015
%   Based on Sam's Aug 4 2015

clear
close all

% N = 50;
N = 120;
epmat = linspace(0.95, 0.999, 10);
epmat = sqrt(0.8:0.01:0.99);

F = load('../../../stokeswavedata/series_steep_N100_ep0.95.mat');
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


solmat = zeros(N+1, length(epmat));

mytol = 1e-10;
options = optimset('Display', 'iter', 'TolFun', mytol, 'TolX', mytol);
for k = 1:length(epmat)
    ep = epmat(k);
    
    if k > 1
        x0 = solmat(:,k-1);
    end
    
    fwd = @(x)vb_series_bern_steep(x, ep, N, usebeta);
    [sol, fval] = fsolve(fwd,x0(:).',options);
    solmat(:,k) = sol;    
end


%% Plot the solution

phi = linspace(-0.5,0.5,5000);
[w, wp] = vb_wfree_steep(phi, solmat(1:N-1, end), solmat(N, end));

zp = 1./w;

x_phi = real(zp);
y_phi = imag(zp);

x = cumtrapz(phi, x_phi) - 0.5;
y = cumtrapz(phi, y_phi);
ycrest = interp1(phi, y, 0);
y = y - ycrest;

figure(1);
plot(x,y); hold all

%% Plot mu vs. ep

figure(2);
plot(epmat, solmat(N+1,:));



