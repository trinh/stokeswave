k=2*pi;
A=0.1;
M=100;
N=5;

phi = linspace(-0.5,0.5,M);

z_an = phi + ((4*1i)./(k*(1+A*exp(-2*pi*1i*phi)))) - (4*1i)/k;

figure(1);
plot(real(z_an), imag(z_an));


est(1) = -8*pi*A;
est(2) = 16*pi*A^2;
est(3) = -24*pi*A^3;
est(4) = 48*pi*A^4;

for j=5:N-1
    est(j) = 0;
end

est(N) = 0.1;
est(N+1) = 1;

z_est = phi;

for j=1:N
    z_est = z_est + ((4*1i)/(k))*((-1)^(j))*(A^j)*exp(-2*pi*1i*j*phi);
end

figure(2);
plot(real(z_est),imag(z_est));

z_dash_est = 1;

for j=1:4
    z_dash_est = z_dash_est + ((-1)^(j))*((8*pi*j)/(k))*(A^j)*(exp(-2*pi*1i.*phi*j));
end

z_dash = 1 - ((8*A*pi)/(k)).*(exp(-2*pi*1i.*phi)).*((1+A.*exp(-2*pi*1i.*phi)).^(-2));

figure(3);
hold on
plot(real(z_dash),imag(z_dash));
plot(real(z_dash_est),imag(z_dash_est));
hold off




