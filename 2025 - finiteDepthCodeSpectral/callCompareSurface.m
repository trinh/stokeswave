%Compare Surfaces
close all
clear

data1 = load("test.mat");
data2 = load("test2.mat");
data3 = load("test3.mat");

Ysam = interp2(real(data1.mesh),imag(data1.mesh),data1.Y,real(data2.ZETA),imag(data2.ZETA));
figure(1)
contourf(real(data2.ZETA),imag(data2.ZETA),Ysam-real(data2.YmatFull)+1)

figure(2)
surf(real(data2.ZETA),imag(data2.ZETA),Ysam,'EdgeColor','none')
hold on
surf(real(data2.ZETA),imag(data2.ZETA),imag(data2.YmatFull),'EdgeColor','none')
%surf(real(data3.ZETA),imag(data3.ZETA),real(data3.YmatFull)-0.5,'EdgeColor','none')