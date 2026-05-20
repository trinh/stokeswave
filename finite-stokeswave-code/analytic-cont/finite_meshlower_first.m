function [Z,Zp] = finite_meshlower_first(wave, F)

phi = real(F(1,:));
psi = imag(F(:,1));
Z = 0*F;
Zp = 0*F;

for j = 1:length(phi)
    current_phi = phi(j);
    for k = 1:length(psi)
        [zz, zpzp] = wave.getZValues(current_phi + 1i*psi(k));
        Z(k, j) = zz; 
        Zp(k, j)= zpzp;
    end
end


%{
figure
% contourf(real(F), imag(F), imag(Z), 50)
p = surf(real(F), imag(F), imag(Z));
set(p, 'EdgeColor', 'none');

figure
% contourf(real(F), imag(F), imag(Z), 50)
p = surf(real(F), imag(F), real(Z));
set(p, 'EdgeColor', 'none');

figure
contourf(real(F), imag(F), real(Z), 50)

figure
contourf(real(F), imag(F), imag(Z), 50)
%}