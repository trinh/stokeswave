% Code to generate analytic continuation paths along arbitrary contours 

addpath(fullfile(pwd,'./paths'))

clear
%close all

% =====================================
%         Load in the solution
F = load('../../stokeswavedata/series_N100_ep0.5.mat');
solseries = F.solseries;
% =====================================


% =====================================
%            Path input
% path = get_path_phil('pub.AA_AA_ep0.5');
path = get_path_2016_08_14('round.0.A');

pathdat = qt_bdint_path_solve(path, solseries);
% pathdat = qt_bdint_path_solve(path, solseries, 'singcoef', [1, 0.17555i]);

omega = pathdat.tau - 1i*pathdat.theta;
W = exp(omega);

%% Plot

% For Henry
% Check the output of the series solution to verify it works with the above
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
% Re-assemble this as omega
omega_ser = log(q_ser) - 1i*theta_ser;
W_ser = exp(omega_ser);

figure(5); clf(5); hold all
plot3(real(pathdat.f),imag(pathdat.f), imag(W), '-', 'LineWidth', 2);
xlabel('Re(f)'); ylabel('Im(f)'); zlabel('Imag W');
plot3(phi, 0*phi, imag(W_ser), '--');
view([ -25.9531   39.9922]);


figure(6); clf(6); hold all
plot3(real(pathdat.f),imag(pathdat.f), real(W), '-', 'LineWidth', 2);
xlabel('Re(f)'); ylabel('Im(f)'); zlabel('Re W');
plot3(phi, 0*phi, real(W_ser), '--');
view([ -25.9531   39.9922]);