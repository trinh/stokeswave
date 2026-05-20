% Code to mesh the surface 
%
% 2026-05-19: Rewrite of the infinite depth version 

clear
close all

% =================================================================
%  Load in the solution 
%  (note data folder should already be in path)
%
%  Solutions are generated via series/getcoef_finite.m script
F = load('series_finitedepth_N101_ep0.7746_r00.3.mat');
wave = F.wave;
% =================================================================

%%

path = get_path_20260520('mesh.0.up');
% For finite depth, need to take care to set PsiMax and PsiMin values
% appropriately since this depends on r0 until we are ready to modify the
% integral code to allow us to step over
pathdat = qt_finite_bdint_mesh_solve(path.f, wave, 'PsiMax', 0.3, 'PsiMin', -0.3);

% path = get_path('mesh_n1_up_test');
% [pathdat_test, ~] = qt_finite_bdint_path_solve(path, solseries, NaN);

%%
%{
path = get_path('mesh_n1_down');
pathdat2 = qt_bdint_mesh_solve(path, solseries);

path = get_path('mesh_n1_down_test');
[pathdat2_test, ~] = qt_bdint_path_solve(path, solseries, NaN);

%%

path = get_path('mesh_n1n1_down');
pathdat3 = qt_bdint_mesh_solve(path, solseries);

path = get_path('mesh_n1n1_down_test');
[pathdat3_test, ~] = qt_bdint_path_solve(path, solseries, NaN);

%%

path = get_path('mesh_n1n1_up');
pathdat4 = qt_bdint_mesh_solve(path, solseries);


%%
%}
%%

% Use 'real', 'imag', 'abs', etc.
reimname = 'imag';
reim = str2func(reimname);

figure(1); clf(1);
hold all
zlabel([reimname, 'tau']);

surf(real(pathdat.FF), imag(pathdat.FF), reim(pathdat.taumat))
plot3(real(pathdat.f), imag(pathdat.f), reim(pathdat.tau), 'k', 'LineWidth', 2)
% plot3(real(pathdat_test.f), imag(pathdat_test.f), reim(pathdat_test.tau), 'b', 'LineWidth', 2)

% Sanity checks 
path_test = [];
path_test.f = [-0.8 + 1e-4*1i, -0.8 + 0.3i, -0.2 + 0.3i, -0.2 + 1e-4*1i];
pathdat_test = qt_finite_bdint_path_solve(path_test.f, wave);
plot3(real(pathdat_test.f), imag(pathdat_test.f), reim(pathdat_test.tau), 'b', 'LineWidth', 2)

% Also plot the wave solution along the axis please
phitest = linspace(-1, 1, 100);
[z, zp] = wave.getZValues(phitest);
[q, theta] = wave.getQTValues(phitest);
plot3(real(phitest), 0*phitest, reim(log(q)), '--', 'LineWidth', 2)


xlabel('phi');
ylabel('psi');
zlabel([reimname, '(tau)']);
title(['ep = ', num2str(wave.ep), ' r0 = ', num2str(wave.r0)]);


%% We can also plot the z values

reimname = 'real';
reim = str2func(reimname);

figure(2); clf(2);
hold all
zlabel([reimname, 'z']);

surf(real(pathdat.FF), imag(pathdat.FF), reim(pathdat.zmat))
surf(real(pathdat2.FF), imag(pathdat2.FF), reim(pathdat2.zmat))
surf(real(pathdat3.FF), imag(pathdat3.FF), reim(pathdat3.zmat))


