% ZPATH_RECTANGLE solves for the z values using the rectangular method of
% analytic continuation

clear
close all


% =====================================
%         Load in the solution
F = load('../../stokeswavedata/series_N100_ep0.6.mat');
solseries = F.solseries;

F = load('../../stokeswavedata/meshdat_M100_ep0.6.mat');
meshdat = F.meshdat;

mu = solseries.mu;
B = solseries.B;
ep = solseries.ep;
% =====================================



mytol = 1e-7;
options = odeset('RelTol', mytol, 'AbsTol', mytol);

% This is the main data structure
recdat = [];

% Initial phi and psi.
phi0 = 0.47;
psi0 = 0.32;

path = [phi0, phi0 - 1i*psi0, -phi0 - 1i*psi0, -phi0, phi0];

M = 50; % 1/4 number of points in the in-fluid mesh. I.e. points on one line
recnum = 4; % No. of rectangles to plot


for rind = 1:recnum 
    recdat(rind).f = [];
    recdat(rind).z = [];
    recdat(rind).zp = [];
    recdat(rind).zminus = [];
    recdat(rind).zpminus = [];
    recdat(rind).P = [];
    for j = 1:length(path)-1
        
        % Path segment end points
        f0 = path(j);
        f1 = path(j+1);
        
        % Compute straight line
        gam = @(s) f0 + s*(f1 - f0);
        dgam = (f1 - f0);
        
        % Use in-fluid data for first one
        if rind == 1
            ss = linspace(0, 1, M);
            ss = ss(:);
            
            recdat(rind).ff{j} = gam(ss);
            [recdat(rind).zz{j}, recdat(rind).zzp{j}] = zfree(gam(ss), solseries);
            recdat(rind).zzminus{j} = [];
            recdat(rind).zzpminus{j} = [];
            recdat(rind).PP{j} = [];
        else
            % Left/right (getting z(-f) and z'(-f) along the previous path segment)
            if j == 1 || j == 3
                Zminus = @(s) recdat(rind-1).Z{j}(-imag(gam(s)));
                Zpminus = @(s) recdat(rind-1).Zp{j}(-imag(gam(s)));
                
                %{
                % Use to test using meshdata for rind == 2
                pU = meshdat.UHP(1).p;
                fU = meshdat.UHP(1).f;
                Zminus = @(s)interp2(real(meshdat.LHP.f), imag(meshdat.LHP.f), meshdat.LHP.z, -real(gam(s)), -imag(gam(s)));
                pfunc = @(s)interp2(real(fU), imag(fU), pU, real(gam(s)), imag(gam(s)));
                Zpminus = pfunc;                
                %}
            elseif j == 2 || j == 4
                Zminus = @(s) recdat(rind-1).Z{j}(-real(gam(s)));
                Zpminus = @(s) recdat(rind-1).Zp{j}(-real(gam(s)));
            end
            
            % Set the initial condition
            if j == 1
                z0 = recdat(rind-1).zz{3}(end);
            else
                z0 = zz(end);
            end
            
            % Solve along the path
            fwd = @(s, z) dgam*zbern_flipode(s, z, solseries, Zminus, Zpminus);
                    
            [ss, zz] = ode45(fwd, [0, 1], z0, options);
            % Re-run to get the P-data
            PP = solseries.mu*solseries.B + 1i*pi*(Zminus(ss) + zz);
            
            zzp = zbern_flipode(ss, zz, solseries, Zminus, Zpminus);
            zzp = zzp(:);
            ss = ss(:);
            zz = zz(:);
            PP = PP(:);
            
            % Save the rectangle data
            recdat(rind).ff{j} = gam(ss);
            recdat(rind).zz{j} = zz;
            recdat(rind).zzp{j} = zzp;
            recdat(rind).zzminus{j} = Zminus(ss);
            recdat(rind).zzpminus{j} = Zpminus(ss);
            recdat(rind).PP{j} = PP;
        end
        
        % Combine the rectangle data along path segments
        recdat(rind).f = [recdat(rind).f; recdat(rind).ff{j}];
        recdat(rind).z = [recdat(rind).z; recdat(rind).zz{j}];
        recdat(rind).zp = [recdat(rind).zp; recdat(rind).zzp{j}];        
        recdat(rind).zminus = [recdat(rind).zminus; recdat(rind).zzminus{j}];
        recdat(rind).zpminus = [recdat(rind).zpminus; recdat(rind).zzpminus{j}];
        recdat(rind).P = [recdat(rind).P; recdat(rind).PP{j}];
    end
    
    % Set up interpolation
    for j = 1:4
        % left/right
        if j == 1 || j == 3
            recdat(rind).Z{j} = @(psi) interp1(imag(recdat(rind).ff{j}), recdat(rind).zz{j}, psi);
            recdat(rind).Zp{j} = @(psi) interp1(imag(recdat(rind).ff{j}), recdat(rind).zzp{j}, psi);
            % up/down
        elseif j == 2 || j == 4
            recdat(rind).Z{j} = @(phi) interp1(real(recdat(rind).ff{j}), recdat(rind).zz{j}, phi);
            recdat(rind).Zp{j} = @(phi) interp1(real(recdat(rind).ff{j}), recdat(rind).zzp{j}, phi);
        end
    end
    
    % Flip the path to repeat.
    path = -path;
    
end

%% Plot
figure(1); figshift

plotme(meshdat.UHP(1).f, real(meshdat.UHP(1).z), 'interp', 10);
hold all
plotme(meshdat.LHP.f, real(meshdat.LHP.z), 'interp', 10);

for j = 1:length(recdat)
    h(j) = plot3(real(recdat(j).f), imag(recdat(j).f), real(recdat(j).z), 'LineWidth', 4);
    lcell{j} = ['j = ', num2str(j)];
end
legend(h, lcell)
xlabel('phi'); ylabel('psi'); zlabel('real(z)');


% figure(2); figshift
% 
% % plotme(meshdat.UHP(1).f, imag(meshdat.UHP(1).z), 'interp', 10); 
% hold all
% % plotme(meshdat.LHP.f, imag(meshdat.LHP.z), 'interp', 10)
% 
% for j = 1:length(recdat)
%     h(j) = plot3(real(recdat(j).f), imag(recdat(j).f), imag(recdat(j).z), 'LineWidth', 3);
%     lcell{j} = ['j = ', num2str(j)];
% end
% 
% xlabel('phi'); ylabel('psi'); zlabel('imag(z)');
% legend(h, lcell)

%%
figure(99); figshift

for j = 2:length(recdat)
    h(j) = plot3(real(recdat(j).f), imag(recdat(j).f), real(recdat(j).P), 'LineWidth', 4); 
    hold all
    lcell{j} = ['j = ', num2str(j)];
end
legend(h, lcell)
xlabel('phi'); ylabel('psi'); zlabel('real(z)');

%%

% figure(3); figshift
% 
% [~,zpLHP] = zfree(meshdat.LHP.f, solseries);
% plotme(meshdat.UHP(1).f, real(meshdat.UHP(1).zp), 'interp', 10); 
% hold all
% plotme(meshdat.LHP.f, real(zpLHP), 'interp', 10);
% 
% for j = 1:length(recdat)
%     h(j) = plot3(real(recdat(j).f), imag(recdat(j).f), real(recdat(j).zp), 'LineWidth', 3);
%     lcell{j} = ['j = ', num2str(j)];
% end
% 
% xlabel('phi'); ylabel('psi'); zlabel('real(zp)');
% legend(h, lcell)

% %%
% figure(4); figshift
% 
% [~,zpLHP] = zfree(meshdat.LHP.f, solseries);
% plotme(meshdat.UHP(1).f, imag(meshdat.UHP(1).zp), 'interp', 10);
% hold all
% plotme(meshdat.LHP.f, imag(zpLHP), 'interp', 10);
% 
% for j = 1:length(recdat)
%     h(j) = plot3(real(recdat(j).f), imag(recdat(j).f), imag(recdat(j).zp), 'LineWidth', 3);
%     lcell{j} = ['j = ', num2str(j)];
% end
% 
% xlabel('phi'); ylabel('psi'); zlabel('real(zp)');
% legend(h, lcell)


%% Plot
% figure(5); figshift; hold all
% for j = 1:length(recdat)
%     h(j) = plot3(real(recdat(j).f), imag(recdat(j).f), abs(1./recdat(j).zp), 'LineWidth', 4);
%     lcell{j} = ['j = ', num2str(j)];
% end
% legend(h, lcell)
% xlabel('phi'); ylabel('psi');

%% For testing purposes, I want to reconstruct the ODE RHS components
% figure(1); figshift; hold all
% for j = 1:length(recdat)
%     if j > 1
%         Q1 = 2*recdat(j).zpminus;
%         Q2 = (solseries.mu*solseries.B + 1i*pi*(recdat(j).z + recdat(j).zminus));
%         RHS = solseries.mu./(Q1.*Q2);
%     end
%     h(j) = plot3(real(recdat(j).f), imag(recdat(j).f), real(recdat(j).zp), 'LineWidth', 4);
%     lcell{j} = ['j = ', num2str(j)];
% end
% legend(h, lcell)
% xlabel('phi'); ylabel('psi');
% title('testing');

%% Save the rectangular data

% name = ['../../stokeswavedata/recdat_rec', num2str(recnum), '_ep', num2str(ep), '.mat'];
% save(name, 'recdat');
