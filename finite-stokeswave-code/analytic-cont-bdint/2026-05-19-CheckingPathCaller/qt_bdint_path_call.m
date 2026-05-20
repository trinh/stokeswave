% Code to generate analytic continuation paths along arbitrary contours 

clear
close all

% =====================================
%         Load in the solution
r0 = 0.5;
ep = sqrt(0.1);
wave = VBStokesWave(r0, ep);
% Plot on a more fine grid
phis =  linspace(-0.5,0.5, 100);
z_ser = [];
zp_ser = [];
for j=1:length(phis)
    [Z, Zp] = wave.getZValues(phis(j));
    zp_ser(j)= Zp;
    z_ser(j) = Z;
end
% q and theta are given by the relationship of log(df/dz) 
Omega = log(1./zp_ser);
tau_ser = real(Omega); q_ser = exp(tau_ser);
theta_ser = -imag(Omega);
% =====================================


% =====================================
%            Path input
path = get_path('round.0.A');
% path = get_path('box1pass');
% path = [];
% path.f = [-0.1 + 1e-4*1i, -0.1 + 0.1i, 0.1 + 0.1i, 0.1 + 1e-4*1i];

% keyboard

pathdat = qt_finite_bdint_path_solve(path, wave);
% pathdat = qt_bdint_path_solve(path, solseries, 'singcoef', [1, 0.17555i]);


%%
omega = pathdat.tau - 1i*pathdat.theta;
W = exp(omega);


figure(5); figshift; hold all
% plot3(real(pathdat.f),imag(pathdat.f), imag(pathdat.omega), '-', 'LineWidth', 2);
% plot3(real(pathdat.f),imag(pathdat.f), imag(W), '-', 'LineWidth', 2);
plot3(real(pathdat.f),imag(pathdat.f), imag(W), '-', 'LineWidth', 2);
xlabel('\phi'); ylabel('\psi');

% Compare with the series solution
plot3(phis, 0*phis, imag(exp(Omega)), '--');

% Run the following code to plot
% handles = qt_bdint_path_plotter(path, pathdat, 'PatchWidth', 0.02, 'NumArc', 150);    