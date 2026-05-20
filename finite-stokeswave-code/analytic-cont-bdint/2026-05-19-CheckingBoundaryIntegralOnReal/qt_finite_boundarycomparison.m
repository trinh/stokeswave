% May 2026
%
% This will be a version similar to the infinite depth
% qt_boundarycomparison that checks the series solutions versus a boundary
% integral computation
%
% There is already a working boundaryintegraltest.m code in the bdint/
% folder but this uses the delta + i*beta formulation instead of the (q,
% theta) formulation

clear
close all

% Load in the series solution and plot
N = 51;         % Should be odd
ep = 0.7746;    % Amplitude
r0 = 0.3;       % Depth parameter

phis =  linspace(-0.5,0.5, N);
wave = VBStokesWave(r0, ep);
z_ser = [];
zp_ser = [];
for j=1:length(phis)
    [Z, Zp] = getZValues(wave, phis(j));
    zp_ser(j)= Zp;
    z_ser(j) = Z;
end

%%
% It is helpful to plot things
figure(1);
plot(real(z_ser), imag(z_ser))
xlabel('x'); ylabel('y');

% q and theta are given by the relationship of log(df/dz) 
Omega = log(1./zp_ser);
tau_ser = real(Omega); q_ser = exp(tau_ser);
theta_ser = -imag(Omega);
% From the series data need to get the speed and angle
% The below is not necessary as the use of Omega is just fine
% q_ser = 1./abs(zp_ser);
% theta_ser = atan(imag(zp_ser)./real(zp_ser));

figure(2);
plot(phis, tau_ser); hold all
plot(phis, theta_ser);
xlabel('\phi'); ylabel('\tau, \theta');

%%

% Now we try and solve on the axis. For inspiration, look into
% qt_boundarycomparison.m from the infinite depth case

bdinit = [q_ser, ...
          theta_ser(2:end-1), ...
          wave.mu, ...
          wave.B];

myeps = 1e-10;
options = optimset('Display', 'iter', 'TolFun', myeps, 'TolX', myeps);
fwd = @(X)qt_finite_bdint_func(X, phis, ep, wave);
[fsol,fval] = fsolve(fwd, bdinit, options); % Call solver
 
N_bd = N;
q_bd = fsol(1:N_bd);
theta_bd(1) = 0;
theta_bd(N_bd) = 0;
for k=2:N_bd-1
    theta_bd(k)=fsol(N_bd+k-1); 
end
mu_bd = fsol(2*N_bd-1);
B_bd = fsol(2*N_bd);

%% Comparison figures

figure(2); 
plot(phis, log(q_bd), '--');
plot(phis, theta_bd, '--');