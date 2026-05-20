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


%%
omega = pathdat.tau - 1i*pathdat.theta;
W = exp(omega);


figure(5); figshift; hold all
% plot3(real(pathdat.f),imag(pathdat.f), imag(pathdat.omega), '-', 'LineWidth', 2);
% plot3(real(pathdat.f),imag(pathdat.f), imag(W), '-', 'LineWidth', 2);
plot3(real(pathdat.f),imag(pathdat.f), imag(W), '-', 'LineWidth', 2);
xlabel('\phi'); ylabel('\psi');

% Run the following code to plot
% handles = qt_bdint_path_plotter(path, pathdat, 'PatchWidth', 0.02, 'NumArc', 150);    