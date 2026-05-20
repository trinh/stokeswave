% Olfe and Rottman method for steep waves
clear
close all

N = 100;
% Initial guess - last one is c

x0(N+1) = sqrt(1.1931); % c

x0(1:N) = zeros(1,N);
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
% Solve the system
options = optimset('Display', 'iter');
fwd = @(x)olfe_eqns(x, N);
[b,~] = fsolve(fwd,x0,options);

%% Plot solution
close all
c = b(N+1);

phi = linspace(-pi,pi,1000);
w = 1;

for k=1:N
    w = w + b(k).*exp(-1i*k.*phi./c);
end

w = (c*(1-exp(-1i.*phi./c)).^(1/3)).*w;
zp = 1./w;


x = cumtrapz(phi,real(zp));
y = cumtrapz(phi,imag(zp));


figure(1);
plot(phi,real(zp));
hold on 
plot(phi,imag(zp));
legend('real','imag');

figure(2);
plot(phi,x);
hold on
plot(phi,y);
legend('real','imag');

figure(3);
plot(x,y);

%{
%% Save all the data
solseries_steep = [];
solseries_steep.N = N;
solseries_steep.an = b(1:N);
solseries_steep.c = c;

name = ['./data/series_olfe_N', num2str(N),  '.mat'];
save(name, 'solseries_steep');
%}
