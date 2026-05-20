% Code to mesh the surface 

clear
close all

% =====================================
%         Load in the solution
F = load('../../stokeswavedata/series_N100_ep0.1.mat');
solseries = F.solseries;
% =====================================

%%

path = get_path('mesh_n1_up');
pathdat = qt_bdint_mesh_solve(path, solseries);

path = get_path('mesh_n1_up_test');
[pathdat_test, ~] = qt_bdint_path_solve(path, solseries, NaN);

%%

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

% Use 'real', 'imag', 'abs', etc.
reimname = 'imag';
reim = str2func(reimname);

figure(1); clf(1);
hold all
zlabel([reimname, 'tau']);

surf(real(pathdat.FF), imag(pathdat.FF), reim(pathdat.taumat))
plot3(real(pathdat.f), imag(pathdat.f), reim(pathdat.tau), 'k', 'LineWidth', 2)
plot3(real(pathdat_test.f), imag(pathdat_test.f), reim(pathdat_test.tau), 'b', 'LineWidth', 2)

xlabel('phi');
ylabel('psi');
zlabel([reimname, '(tau)']);
title(['ep = ', num2str(solseries.ep)]);

%%
surf(real(pathdat2.FF), imag(pathdat2.FF), reim(pathdat2.taumat))
plot3(real(pathdat2.f), imag(pathdat2.f), reim(pathdat2.tau), 'k', 'LineWidth', 2)
plot3(real(pathdat2_test.f), imag(pathdat2_test.f), reim(pathdat2_test.tau), 'g', 'LineWidth', 2)

%%

surf(real(pathdat3.FF), imag(pathdat3.FF), reim(pathdat3.taumat))
plot3(real(pathdat3.f), imag(pathdat3.f), reim(pathdat3.tau), 'k', 'LineWidth', 2)
plot3(real(pathdat3_test.f), imag(pathdat3_test.f), reim(pathdat3_test.tau), 'm', 'LineWidth', 2)


%% We can also plot the z values

reimname = 'real';
reim = str2func(reimname);

figure(2); clf(2);
hold all
zlabel([reimname, 'z']);

surf(real(pathdat.FF), imag(pathdat.FF), reim(pathdat.zmat))
surf(real(pathdat2.FF), imag(pathdat2.FF), reim(pathdat2.zmat))
surf(real(pathdat3.FF), imag(pathdat3.FF), reim(pathdat3.zmat))


