% Code to mesh the surface
% like qt_meshsurface.m but for individual paths
%
% ******** Generates for epsilon = 0.3

clear
close all

% =====================================
%         Load in the solution
F = load('../../stokeswavedata/series_N100_ep0.3.mat');
solseries = F.solseries;
% =====================================

exportme = 'true';

psimax = 1.5;
psimin = -1.5;
npsi = 50;
nphi = 100;

pathlist = {'dualmesh.0', ... % 1
            'dualmesh.A', ... % 2
            'dualmesh.AA_', ... % 3
            'dualmesh.AA_B_', ...
            'dualmesh.AA_B_^', ...
            'dualmesh.A1', ...
            'dualmesh.An1', ...
            'dualmesh.AA_AA_', ...
            'dualmesh.AA_AA_B_', ...
            'dualmesh.AA_AA_B_^', ...
            'dualmesh.AA_AA_C_', ...
            'dualmesh.AA_AA_C_^', ...
            'dualmesh.A1A1_', ...
            'dualmesh.A1A_', ...
            'dualmesh.AA_B_^B_'};
        
pathlist = pathlist(end);

for j = 1:length(pathlist)
    pathname = pathlist{j};
    name = ['../../stokeswavepub/contourep0p3/', pathname, '.mat'];    
    
    %% Load the paths
    path = get_path_ep0p3(pathname);
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
    
    drawnow
    
    switch exportme
        case 'true'
            save(name);
    end
end