% Get the coefficients of z'(f) for steep waves
% -----------------------------
% Written July 26 2015
% From Vanden-Broeck, 1986 Phys. Fluids

clear
close all

N = 100;
ep = sqrt(0.6);
%{
F = load('series_steep_N100_ep0.7746_usebeta0.mat');
b = F.steepsolseries.bn;
mu = F.steepsolseries.mu;
% b = convert_z2w(solseries);
%}

usebeta = 1; % 1 to use beta formulation
%initial guess
x0 = zeros(1,N);
if usebeta == 0
    x0(1:N) = b(1:N);
    x0(N+1) = mu;
elseif usebeta == 1    
    x0(1:N-1) = b(1:N-1);
    x0(N) = 1e-5; % This is beta
    x0(N+1) = mu;   
end

% F = load('series_steep_N60_ep0.1.mat');
% x0(1:N-1) = F.steepsolseries.bn;
% x0(N) = F.steepsolseries.mu;
% x0(N+1) = F.steepsolseries.beta;

mytol = 1e-6;
options = optimset('Display', 'iter', 'TolFun', mytol, 'TolX', mytol);
fwd = @(x)series_bern_steep2(x, ep, N, usebeta);
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
% keyboard

% name = ['series_steep_N', num2str(N), '_ep', num2str(ep), '_usebeta', num2str(usebeta), '.mat'];
% save(name, 'steepsolseries');



%% Plot the solution

phi = linspace(-0.5,0.5,1000);
[w, wp] = wfree_steep(phi, steepsolseries.bn, steepsolseries.beta);
% w = 1;
% for k=1:N-1
%     w = w + est(k).*(exp(-2*pi*1i*phi));
% end
% 
% w = ((1-est(N+1).*exp(-2*pi*1i.*phi)).^(1/3)).*w;

zp = 1./w;

x_phi = real(zp);
y_phi = imag(zp);

x = cumtrapz(phi, x_phi) - 0.5;
y = cumtrapz(phi, y_phi);
ycrest = interp1(phi, y, 0);
y = y - ycrest;

figure(1);
plot(x,y); hold all
% zt = zfree(phi, solseries);
% plot(real(zt), imag(zt), 'bo');



