clear

% =====================================
%         Load in the solution
F = load('../data/series_N100_ep0.05.mat');
% F = load('data/series_N100_ep0.7746.mat');
% F = load('data/series_N101_ep0.7.mat');
solseries = F.solseries;
N = solseries.N;
ep = solseries.ep;
% =====================================

% x0 = zeros(1,2*N-2);
% x0(2*N-1) = 1;
% x0(2*N) = 1;
% 
% delta_guess_interp = interp1(linspace(-0.5,0.5,51), delta_guess, linspace(-0.5,0.5,N));
% beta_guess_interp = interp1(linspace(-0.5,0.5,51), beta_guess, linspace(-0.5,0.5,N));
% x0(1:N) = delta_guess_interp;
% x0(N+1:2*N-2) = beta_guess_interp(2:N-1);

N_bd = 201;
phi = linspace(-0.5, 0.5, N_bd);

% Get series z
[z_ser, zp_ser] = zfree(phi, solseries);
z_ser = z_ser(:).'; zp_ser = zp_ser(:).';

bdinit = [real(zp_ser), ...
          imag(zp_ser(2:end-1)), ...
          solseries.mu, ...
          solseries.B];

myeps = 1e-10;
options = optimset('Display', 'iter', 'TolFun', myeps, 'TolX', myeps);
fwd = @(X)boundaryintegral(X, phi, ep);
[fsol,fval] = fsolve(fwd, bdinit, options); % Call solver

delta_bd = fsol(1:N_bd);
beta_bd(1) = 0;
beta_bd(N_bd) = 0;
for k=2:N_bd-1
    beta_bd(k)=fsol(N_bd+k-1); 
end
mu_bd = fsol(2*N_bd-1);
B_bd = fsol(2*N_bd);

% Get X and Y
m = (N_bd + 1)/2;
dphi = phi(2) - phi(1);
Y_bd = 0*delta_bd;
for j = m:N_bd - 1
    Y_bd(j+1) = Y_bd(j) + dphi/2*(beta_bd(j+1) + beta_bd(j));
end
for j = m:-1:2
    Y_bd(j-1) = Y_bd(j) - dphi/2*(beta_bd(j-1) + beta_bd(j));
end
X_bd = cumtrapz(phi, delta_bd) - 0.5;

solbdint = [];
solbdint.delta = delta_bd;
solbdint.beta = beta_bd;
solbdint.mu = mu_bd;
solbdint.B = B_bd;
solbdint.x = X_bd;
solbdint.y = Y_bd;
solbdint.phi = phi;
solbdint.ep = ep;
solbdint.N = N_bd;
% name = ['./data/bdint_N', num2str(N_bd), '_ep', num2str(ep), '.mat'];
% save(name, 'solbdint');

%%

% figure(1);
% plot(phi, delta_bd); hold all
% plot(phi, real(zp_ser), '--');

figure(2);
plot(X_bd, Y_bd); hold all
plot(real(z_ser), imag(z_ser), '--');
xlabel('x'); ylabel('y');

%%

bernerror = 1/2./(delta_bd.^2 + beta_bd.^2) + 2*pi/solseries.mu*Y_bd - B_bd;
% figure(3);
% plot(phi, bernerror);