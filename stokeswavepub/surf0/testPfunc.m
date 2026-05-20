clear
close all

% Testing whether P function is zero
F = load('mesh0_ep0.3_100100_sym.mat');
G = load('../../stokeswavedata/series_N100_ep0.3.mat');
meshdat = F.pathdat;

fup = meshdat.up.FF;
zup = meshdat.up.zmat;
tauup = meshdat.up.taumat;
thetaup = meshdat.up.thetamat;
omegaup = tauup - 1i*thetaup;

fdown = meshdat.down.FF;
zdown = meshdat.down.zmat;
taudown = meshdat.down.taumat;
thetadown = meshdat.down.thetamat;
omegadown = taudown - 1i*thetadown;

B = G.solseries.B;
mu = G.solseries.mu;

P = mu*B + pi*1i*zdown;
Pflip = P;

[C, h] = contourf(real(fup), imag(fup), abs(Pflip + 1i*pi*zup), 50);

% %%
% figure
% [C, h] = contourf(real(fdown), imag(fdown), abs(exp(omegadown)), 50);
% hold on
% [C, h] = contourf(real(fup), imag(fup), abs(exp(omegaup)), 50);
%%
figure
% [C, h] = contourf(real(fdown), imag(fdown), real(zdown), 50);
% hold on
[C, h] = contourf(real(fup), imag(fup), real(zup), 50);
%%
% figure
% [C, h] = contourf(real(fdown), imag(fdown), abs(1./exp(3*taudown)), 50);