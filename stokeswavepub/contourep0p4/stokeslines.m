clear

load dualmesh.0.mat;

w = pathdat.up.FF;
phi = real(w(1,:));
psi = imag(w(:,1));
tau = pathdat.up.taumat;
theta = pathdat.up.thetamat;

q0 = exp(tau);

A = 0.131i;

myint = 0*q0;
myint(1,:) = cumtrapz(phi, q0(1,:));
for k = 1:size(w,2)
    tmp = 1i*cumtrapz(psi, q0(:,k));
    myint(:,k) = myint(1,k) + tmp;
end
chi = -1i*myint;
chic = interp2(real(w), imag(w), chi, real(A), imag(A));
chi = chi - chic;

% [C, h] = contourf(real(w), imag(w), imag(chi), 50);
[C, h] = contourf(real(w), imag(w), real(chi), 50);
hold on
[C, h] = contour(real(w), imag(w), imag(chi), [0 0],'m');
hold off
colorbar