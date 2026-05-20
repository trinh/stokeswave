clear 
close all

r0 = 0;     %Depth
ep = 0.9;   %Steepness
N = 101;    %Number of Sample points

%construct wave
phis =  linspace(-0.5,0.5, N);
    wave = VBStokesWave(r0, ep);
    z = [];
    zp = [];
        for j=1:length(phis)
            [Z, Zp] = getZValues(wave, phis(j));
            zp(j)= Zp;
            z(j) = Z;
        end


figure(1)
plot(real(z),imag(z));
xlabel('real(z)'); ylabel('imag(z)');

figure(2); figshift
plot(real(zp),imag(zp));
xlabel('real(zp)'); ylabel('imag(zp)');

figure(3); figshift
plot(phis,abs(z),'-'); hold all
plot(phis,abs(zp), '--');
xlabel('phi'); 
legend('abs(z)', 'abs(zp)');

figure(4); figshift
plot(phis,real(zp),'-'); hold all
plot(phis,imag(zp), '--');
xlabel('phi'); 
legend('real(zp)', 'imag(zp)');
