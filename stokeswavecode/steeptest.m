clear
close all

F = load('data/series_steep_N100_ep0.99.mat');
steepsolseries = F.steepsolseries;
bn = steepsolseries.bn;
beta = steepsolseries.beta;

L = 50;
M=50;

phi = linspace(-0.5,0.5,M);
psi = linspace(-0.5,0,M);
[Phi,Psi] = meshgrid(phi,psi);
f = Phi+1i*Psi;

for k1=1:M
for k2=1:M
[z(k1,k2),zp(k1,k2)] = zfree_steep(f(k1,k2), steepsolseries, L);
end
end

figure(1);
surf(real(f),imag(f),real(zp));
figure(2);
surf(real(f),imag(f),imag(zp));
figure(3);
surf(real(f),imag(f),real(z));
figure(4);
surf(real(f),imag(f),imag(z));
