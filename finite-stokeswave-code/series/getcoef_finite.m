% Generate the series solutions for finite depth and save them in an
% appropriate directory 
% 
% Written May 2026; modelled after getcoef in the infinite depth 

clear
close all

% Load in the series solution and plot
N = 101;         % Should be odd
ep = 0.5;    % Amplitude
r0 = 0.01;       % Depth parameter

wave = VBStokesWave(r0, ep);

name = ['./../../finite-stokeswave-data/series_finitedepth_N', num2str(N), ...
        '_ep', num2str(ep), ...
        '_r0', num2str(r0), ...
        '.mat'];
save(name, 'wave');

%%

phis =  linspace(-0.5,0.5, 200);
z_ser = [];
zp_ser = [];
for j=1:length(phis)
    [Z, Zp] = getZValues(wave, phis(j));
    zp_ser(j)= Zp;
    z_ser(j) = Z;
end


% It is helpful to plot things
figure(1);
plot(real(z_ser), imag(z_ser))
xlabel('x'); ylabel('y');

% q and theta are given by the relationship of log(df/dz) 
Omega = log(1./zp_ser);
tau_ser = real(Omega); q_ser = exp(tau_ser);
theta_ser = -imag(Omega);

