clear
close all

F = load('../../stokeswavedata/series_steep_N100_ep0.8.mat');
steepsolseries = F.steepsolseries;


phi = linspace(-0.5,0.5,25);
psi = linspace(-0.5,0,25);
[Phi,Psi] = meshgrid(phi,psi);
L = 50;



z=[];

for k1=1:25
for k2=1:25
    z(k1,k2) = zfree_steep(Phi(k1,k2) +1i*Psi(k1,k2), steepsolseries, L);
end
end

figure(1);
surf(Phi,Psi,imag(z));