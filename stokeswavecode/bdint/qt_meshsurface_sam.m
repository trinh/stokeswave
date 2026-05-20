% Code to mesh the surface 

clear
close all

% =====================================
%         Load in the solution
F = load('../../stokeswavedata/series_N100_ep0.6.mat');
solseries = F.solseries;
% =====================================

% May not want to run all these at once
% Need to update the notation

%{
path = get_path('mesh_0_up');
pathdat = qt_bdint_mesh_solve(path, solseries);

path = get_path('mesh_B_up');
pathdat2 = qt_bdint_mesh_solve(path, solseries);

path = get_path('mesh_B_down');
pathdat3 = qt_bdint_mesh_solve(path, solseries);

path = get_path('mesh_BE_down');
pathdat4 = qt_bdint_mesh_solve(path, solseries);

path = get_path('mesh_BE_up');
pathdat5 = qt_bdint_mesh_solve(path, solseries);

path = get_path('mesh_BED_down');
pathdat6 = qt_bdint_mesh_solve(path, solseries);

path = get_path('mesh_BED_up');
pathdat7 = qt_bdint_mesh_solve(path, solseries);

path = get_path('mesh_BEDB_up');
pathdat8 = qt_bdint_mesh_solve(path, solseries);

path = get_path('mesh_BEDB_down');
pathdat9 = qt_bdint_mesh_solve(path, solseries);
%}

path = get_path('mesh_BEBE_down');
pathdat_tmp = qt_bdint_mesh_solve(path, solseries);

%{
% Example 3
path = get_path('mesh_BEDE_down');
pathdat_extra_down = qt_bdint_mesh_solve(path, solseries);

path = get_path('mesh_BEDE_up');
pathdat_extra_up = qt_bdint_mesh_solve(path, solseries);
%}

% Figure examples
%{
path = get_path('mesh_aAB_down');
pathdat_aAB_down = qt_bdint_mesh_solve(path, solseries);
%}

%{
path = get_path('mesh_aABlab_up');
pathdat_aABlab_up = qt_bdint_mesh_solve(path, solseries);

path = get_path('mesh_aABlab_down');
pathdat_aABlab_down = qt_bdint_mesh_solve(path, solseries);

path = get_path('mesh_aAaABlA_down');
pathdat_aAaABlA_down = qt_bdint_mesh_solve(path, solseries);

%}

%%
% Use 'real', 'imag', 'abs', etc.
reimname = 'imag';
reim = str2func(reimname);


% Each figure represents one sheet

%{
figure(1); hold all 
surf(real(pathdat.FF), imag(pathdat.FF), reim(pathdat.taumat))
plot3(real(pathdat.f), imag(pathdat.f), reim(pathdat.tau), 'k', 'LineWidth', 2)

figure(2); hold all 
surf(real(pathdat2.FF), imag(pathdat2.FF), reim(pathdat2.taumat))
plot3(real(pathdat2.f), imag(pathdat2.f), reim(pathdat2.tau), 'k', 'LineWidth', 2)

surf(real(pathdat3.FF), imag(pathdat3.FF), reim(pathdat3.taumat))
plot3(real(pathdat3.f), imag(pathdat3.f), reim(pathdat3.tau), 'k', 'LineWidth', 2)

figure(3); hold all 
surf(real(pathdat4.FF), imag(pathdat4.FF), reim(pathdat4.taumat))
plot3(real(pathdat4.f), imag(pathdat4.f), reim(pathdat4.tau), 'k', 'LineWidth', 2)

surf(real(pathdat5.FF), imag(pathdat5.FF), reim(pathdat5.taumat))
plot3(real(pathdat5.f), imag(pathdat5.f), reim(pathdat5.tau), 'k', 'LineWidth', 2)

figure(4); hold all 
surf(real(pathdat6.FF), imag(pathdat6.FF), reim(pathdat6.taumat))
plot3(real(pathdat6.f), imag(pathdat6.f), reim(pathdat6.tau), 'k', 'LineWidth', 2)

surf(real(pathdat7.FF), imag(pathdat7.FF), reim(pathdat7.taumat))
plot3(real(pathdat7.f), imag(pathdat7.f), reim(pathdat7.tau), 'k', 'LineWidth', 2)

figure(7); hold all
surf(real(pathdat8.FF), imag(pathdat8.FF), reim(pathdat8.taumat))
plot3(real(pathdat8.f), imag(pathdat8.f), reim(pathdat8.tau), 'k', 'LineWidth', 2)

surf(real(pathdat9.FF), imag(pathdat9.FF), reim(pathdat9.taumat))
plot3(real(pathdat9.f), imag(pathdat9.f), reim(pathdat9.tau), 'k', 'LineWidth', 2)

%  BEDE example
figure(6); hold all
surf(real(pathdat_extra_down.FF), imag(pathdat_extra_down.FF), reim(pathdat_extra_down.taumat))
plot3(real(pathdat_extra_down.f), imag(pathdat_extra_down.f), reim(pathdat_extra_down.tau), 'k', 'LineWidth', 2)

surf(real(pathdat_extra_up.FF), imag(pathdat_extra_up.FF), reim(pathdat_extra_up.taumat))
plot3(real(pathdat_extra_up.f), imag(pathdat_extra_up.f), reim(pathdat_extra_up.tau), 'k', 'LineWidth', 2)
%}

% Tmp test
figure(7); hold all
surf(real(pathdat_tmp.FF), imag(pathdat_tmp.FF), reim(pathdat_tmp.taumat))
plot3(real(pathdat_tmp.f), imag(pathdat_tmp.f), reim(pathdat_tmp.tau), 'k', 'LineWidth', 2)


%{
figure(8); hold all
surf(real(pathdat_aAB_down.FF), imag(pathdat_aAB_down.FF), reim(pathdat_aAB_down.taumat))
plot3(real(pathdat_aAB_down.f), imag(pathdat_aAB_down.f), reim(pathdat_aAB_down.tau), 'k', 'LineWidth', 2)


figure(9); hold all
surf(real(pathdat_aABlab_up.FF), imag(pathdat_aABlab_up.FF), reim(pathdat_aABlab_up.taumat))
plot3(real(pathdat_aABlab_up.f), imag(pathdat_aABlab_up.f), reim(pathdat_aABlab_up.tau), 'k', 'LineWidth', 2)
surf(real(pathdat_aABlab_down.FF), imag(pathdat_aABlab_down.FF), reim(pathdat_aABlab_down.taumat))
plot3(real(pathdat_aABlab_down.f), imag(pathdat_aABlab_down.f), reim(pathdat_aABlab_down.tau), 'k', 'LineWidth', 2)


figure(10); hold all
surf(real(pathdat_aAaABlA_down.FF), imag(pathdat_aAaABlA_down.FF), reim(pathdat_aAaABlA_down.taumat))
plot3(real(pathdat_aAaABlA_down.f), imag(pathdat_aAaABlA_down.f), reim(pathdat_aAaABlA_down.tau), 'k', 'LineWidth', 2)
%}


xlabel('phi');
ylabel('psi');
zlabel([reimname, '(tau)']);
title(['ep = ', num2str(solseries.ep)]);


