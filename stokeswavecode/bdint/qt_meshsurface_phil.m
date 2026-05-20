% Code to mesh the surface 
% like qt_meshsurface.m but for individual paths

clear
close all

% =====================================
%         Load in the solution
F = load('../../stokeswavedata/series_N100_ep0.5.mat');
solseries = F.solseries;
% =====================================

psimax = 0.6;
psimin = -0.8;
npsi = 100;
nphi = 200;

%% Load the paths
path = get_path_ep0p5('dualmesh.AA_');
pathdat.up = qt_bdint_mesh_solve(path.up, solseries, 'PsiMax', psimax, 'PsiMin', psimin, 'NPsi', npsi, 'NPhi', nphi);
pathdat.down = qt_bdint_mesh_solve(path.down, solseries, 'PsiMax', psimax, 'PsiMin', psimin, 'NPsi', npsi, 'NPhi', nphi);

%% Run an extra path
pathdat.test = qt_bdint_path_solve(path, solseries);


%%
% Use 'real', 'imag', 'abs', etc.
reimname = 'imag';
reim = str2func(reimname);

figure(1); clf(1);
hold all
zlabel([reimname, 'tau']);

surf(real(pathdat.down.FF), imag(pathdat.down.FF), reim(exp(pathdat.down.taumat - 1i*pathdat.down.thetamat)))
surf(real(pathdat.up.FF), imag(pathdat.up.FF), reim(exp(pathdat.up.taumat - 1i*pathdat.up.thetamat)))
plot3(real(pathdat.test.f), imag(pathdat.test.f), reim(exp(pathdat.test.omega)), 'b', 'LineWidth', 2)

xlabel('phi');
ylabel('psi');
zlabel([reimname, '(tau)']);
title(['ep = ', num2str(solseries.ep)]);


%% We can also plot the z values
%{
reimname = 'real';
reim = str2func(reimname);

figure(2); clf(2);
hold all
title([reimname, '(z)']);

surf(real(pathdat.down.FF), imag(pathdat.down.FF), reim(pathdat.down.zmat))
surf(real(pathdat.up.FF), imag(pathdat.up.FF), reim(pathdat.up.zmat))
%}

%% We can also plot the zp values (note df/dz = u - iv = exp(omega))
%{
reimname = 'real';
reim = str2func(reimname);

figure(3); clf(3);
hold all
title([reimname, '(dz/df)']);

omegamat = pathdat.down.taumat - 1i*pathdat.down.thetamat;
surf(real(pathdat.down.FF), imag(pathdat.down.FF), reim(exp(-omegamat)))
omegamat = pathdat.up.taumat - 1i*pathdat.up.thetamat;
surf(real(pathdat.up.FF), imag(pathdat.up.FF), reim(exp(-omegamat)))
%}

%% Convert to polar

% L = 1;
% seg = ceil((xmax-L/2)/L);
% 
% for j = -seg:seg
%     X = real(pathdat.down.FF);
%     Y = imag(pathdat.down.FF);
%     Z = reim(pathdat.down.taumat);
%     xlow = j*L - L/2;
%     xhigh = xlow + 1;
%     Y((X < xlow) | (X > xhigh)) = NaN;
%     Z((X < xlow) | (X > xhigh)) = NaN;
%     T = exp(-2i*pi*(X + 1i*Y));
%     figure(12); clf(12);
%     contourf(real(T), imag(T), Z, 50);
%     keyboard
% end
% 
% % We need to first split into segments
% xmax = real(pathdat.up.FF(1,end));
% xmin = real(pathdat.up.FF(1,1));

