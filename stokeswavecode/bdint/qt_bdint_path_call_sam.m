% Code to generate analytic continuation paths along arbitrary contours 
% _sam version does singularity location and coefficient

clear
%close all


% ====================================
%         Load in the solution
F = load('../../stokeswavedata/series_N100_ep0.4.mat');
solseries = F.solseries;
% =====================================


% =====================================
%            Path input
path = get_path('b-aAB^a');

sing_start=13; % The point at which the path looping around the singularity starts. 

[~, singloc] = qt_bdint_path_solve_sam_singloc(path, solseries, sing_start); % Gets the singularity location
[pathdat, singcoeff] = qt_bdint_path_solve_sam_singcoeff(path, solseries, sing_start, -singloc); % Gets the path and the singularity coefficient

%%
%W = exp(pathdat.omega);
W = exp(pathdat.tau - 1i*pathdat.theta);



figure(4); figshift; hold all
plot3(real(pathdat.f),imag(pathdat.f), imag(1./W), '-', 'LineWidth', 2);
% plot3(real(pathdat.f),imag(pathdat.f), imag(pathdat.theta), 'b--', 'LineWidth', 2);
xlabel('\phi'); ylabel('\psi');


%% Singularity location

-singloc
singcoeff
