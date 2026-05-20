clear
close all

load('series_N200_ep0.6.mat')


phimin = -1;
phimax = 1;
psimin = -0.5;
psimax = 0;

% Generate N particles
N = 4000;
phi_p = phimin + (phimax - phimin)*rand(N, 1);
psi_p = psimin + (psimax - psimin)*rand(N, 1);

% psi_p = linspace(psimin, psimax, 500);
% phi_p = 0*psi_p;

f = phi_p + 1i*psi_p;
f = f(:);

f_l = linspace(phimin, phimax, 100);
f_l = f_l(:);

[z, zp] = flowfield_zfree(f, solseries);
[z_l, zp_l] = flowfield_zfree(f_l, solseries);
q = 1./zp;
q_l = 1./zp_l;

figure
hold on
circsize = 10;
p = scatter(real(z), imag(z), circsize, imag(z), 'filled');
pl = plot(real(z_l), imag(z_l), 'k.-');
% xlim([-0.5, 0.5]);

dt = 1e-2;
c = .5;
for j = 1:500
    f = f + dt*(abs(q).^2);
    f_l = f_l + dt*(abs(q_l).^2);
    
%     ind = find(real(f) > phimax);
%     newpsi = psimin + (psimax - psimin)*rand(length(ind), 1);
%     f(ind) = phimin + 1i*newpsi(:);
%     
%     ind = find(real(f_l) > phimax);
%     f_l(ind) = f_l(ind) - (phimax - phimin);
%     f_l = sort(f_l);
    
    [z, zp] = flowfield_zfree(f, solseries);
    [z_l, zp_l] = flowfield_zfree(f_l, solseries);
    
    q = 1./zp;
    q_l = 1./zp_l;
%     p = scatter(real(z)-j*dt, imag(z), circsize, imag(z), 'filled');
    set(p, 'XData', real(z) - c*j*dt, 'YData', imag(z), 'CData', imag(z));
    set(pl, 'XData', real(z_l) - c*j*dt, 'YData', imag(z_l));
    drawnow;
end