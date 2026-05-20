k = 6;
A = 0.4;
M = 100;

phi = linspace(0,1,M);
psi = linspace(-0.5,0.5,M);

[Phi, Psi] = meshgrid(phi,psi);

f = Phi+1i*Psi;

z = f + ((4*1i)./(k*(1+A*exp(-2*1i*pi*f)))) - ((4*1i)/(k));

figure(1);
surf(real(z),imag(z),Psi);

figure(2);
surf(Phi,Psi,real(z));

figure(3);
surf(Phi,Psi,imag(z));

figure(4);
contourf(Phi,Psi,imag(z));

phi_1d = linspace(-0.5,0.5);

z_1d = phi_1d + ((4*1i)./(k*(1+A*exp(-2*1i*pi*phi_1d)))) - ((4*1i)/(k));

figure(5);
plot(real(z_1d),imag(z_1d));
axis equal


syms f1
z2 = symfun(f1 + ((4*1i)./(k*(1+A*exp(-2*1i*pi*f1)))) - ((4*1i)/(k)),f1);



