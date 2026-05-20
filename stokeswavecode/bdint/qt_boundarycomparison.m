% QT_BOUNDARYCOMPARISON
% Implements the boundary integral solver comparing to series data
% (q, theta) formulation using the Bernoulli equation with B unknown

% 18 May 2026
% Re-edit to see if Henry can be shown correctedness. Added the flag to not
% save to file
%
% The purpose of this script is to check that the series solution agrees
% with the boundary integral solution on the real axis

clear
savefile = 0; 

% =====================================
%         Load in the solution
F = load('./../../stokeswavedata/series_N100_ep0.2.mat');
solseries = F.solseries;
N = solseries.N;
ep = solseries.ep;
% =====================================

N_bd = 201;
phi = linspace(-0.5, 0.5, N_bd);

% Get series z
[z_ser, zp_ser] = zfree(phi, solseries);
z_ser = z_ser(:).'; zp_ser = zp_ser(:).';

% From the series data need to get the speed and angle
q_ser = 1./abs(zp_ser);
theta_ser = atan(imag(zp_ser)./real(zp_ser));

bdinit = [q_ser, ...
          theta_ser(2:end-1), ...
          solseries.mu, ...
          solseries.B];

myeps = 1e-10;
options = optimset('Display', 'iter', 'TolFun', myeps, 'TolX', myeps);
fwd = @(X)qt_bdint_func(X, phi, ep);
[fsol,fval] = fsolve(fwd, bdinit, options); % Call solver

q_bd = fsol(1:N_bd);
theta_bd(1) = 0;
theta_bd(N_bd) = 0;
for k=2:N_bd-1
    theta_bd(k)=fsol(N_bd+k-1); 
end
mu_bd = fsol(2*N_bd-1);
B_bd = fsol(2*N_bd);

% Get X and Y
[X_bd, Y_bd] = qt_getxy(phi, q_bd, theta_bd);

solqtbdint = [];
solqtbdint.q = q_bd;
solqtbdint.theta = theta_bd;
solqtbdint.mu = mu_bd;
solqtbdint.B = B_bd;
solqtbdint.x = X_bd;
solqtbdint.y = Y_bd;
solqtbdint.phi = phi;
solqtbdint.ep = ep;
solqtbdint.N = N_bd;

% Added May 2026; note this will break; I don't want to overwrite data for
% now
if savefile == 1
    name = ['./../data/bdint_qt_N', num2str(N_bd), '_ep', num2str(ep), '.mat'];
    save(name, 'solqtbdint');
end

%%

figure(2);
plot(X_bd, Y_bd); hold all
plot(real(z_ser), imag(z_ser), '--');
xlabel('x'); ylabel('y');

figure(3);
plot(phi, theta_bd); hold all
plot(phi, theta_ser, 'o');
