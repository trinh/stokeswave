% Compare boundary integral in-fluid mesh with zfree
clear
close all

% =====================================
%         Load in the solution
F = load('./../../stokeswavedata/bdint_N801_ep0.05.mat');
solbdint = F.solbdint;
F = load('./../../stokeswavedata/series_N100_ep0.05.mat');
solseries = F.solseries;
F = load('./../../stokeswavedata/meshdatbd_ep0.05.mat');
meshdatbd = F.meshdatbd;

wmat = meshdatbd.wmat;
[z_ser, zp_ser] = zfree(meshdatbd.wmat, solseries);

z_bd = meshdatbd.xmat + 1i*meshdatbd.ymat;
zp_bd = meshdatbd.xpmat + 1i*meshdatbd.ypmat;
surf(real(wmat), imag(wmat), abs(zp_ser - zp_bd));
% surf(real(wmat), imag(wmat), abs(z_ser - z_bd));