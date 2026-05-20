clear
close all


% =====================================
%         Load in the solution
F = load('../../stokeswavedata/series_steep_N100_ep0.995.mat');
steepsolseries = F.steepsolseries;

mu = steepsolseries.mu;
ep = steepsolseries.ep;
% =====================================

mytol = 1e-11;
options = odeset('RelTol', mytol, 'AbsTol', mytol);

% This is the main data structure
recdat = [];

% Initial phi and psi.
phi0 = 0.1;
psi0 = 0.0021;

path = [phi0, phi0 - 1i*psi0, -phi0 - 1i*psi0, -phi0];
path = [path path phi0];

M = 50; % 1/4 number of points in the in-fluid mesh. I.e. points on one line
recnum = 2; % No. of rectangles to plot
L = 1000;

for rind = 1:recnum 
    recdat(rind).f = [];
    recdat(rind).z = [];
    recdat(rind).zp = [];
    recdat(rind).g = [];
    recdat(rind).gp = [];
    
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
            [recdat(rind).zz{j}, recdat(rind).zzp{j}] = zfree_steep(gam(ss),steepsolseries, L);
        else
            % Left/right (getting z(-f) and z'(-f) along the previous path segment)
            if mod(j,4) == 1 || mod(j,4) == 3
                Zminus = @(s) recdat(rind-1).Z{j}(-imag(gam(s)));
                Zpminus = @(s) recdat(rind-1).Zp{j}(-imag(gam(s)));
            % Up/down
            elseif mod(j,4) == 2 || mod(j,4) == 0
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
            fwd = @(s, z) dgam*zbern_flipode_steep(s, z, steepsolseries, Zminus, Zpminus);
                    
%             [ss, zz] = ode45(fwd, [0, 1], z0, options);
            [ss, zz] = ode113(fwd, [0, 1], z0, options);
            
            zzp = zbern_flipode_steep(ss, zz, steepsolseries, Zminus, Zpminus);
            zzp = zzp(:);
            ss = ss(:);
            zz = zz(:);
            
            % Save the rectangle data
            recdat(rind).ff{j} = gam(ss);
            recdat(rind).zz{j} = zz;
            recdat(rind).zzp{j} = zzp;
            recdat(rind).gg{j} = zz + Zminus(ss);
            recdat(rind).ggp{j} = zzp - Zpminus(ss);
        end
        
        % Combine the rectangle data along path segments
        recdat(rind).f = [recdat(rind).f; recdat(rind).ff{j}];
        recdat(rind).z = [recdat(rind).z; recdat(rind).zz{j}];
        recdat(rind).zp = [recdat(rind).zp; recdat(rind).zzp{j}];
        if rind == 1
        else
        recdat(rind).g = [recdat(rind).g; recdat(rind).gg{j}];
        recdat(rind).gp = [recdat(rind).gp; recdat(rind).ggp{j}];   
        end
    end
    
    % Set up interpolation
    for j = 1:length(path)-1
        % left/right
        if mod(j,4) == 1 || mod(j,4) == 3
            recdat(rind).Z{j} = @(psi) interp1(imag(recdat(rind).ff{j}), recdat(rind).zz{j}, psi);
            recdat(rind).Zp{j} = @(psi) interp1(imag(recdat(rind).ff{j}), recdat(rind).zzp{j}, psi);
        % up/down
        elseif mod(j,4) == 2 || mod(j,4) == 0
            recdat(rind).Z{j} = @(phi) interp1(real(recdat(rind).ff{j}), recdat(rind).zz{j}, phi);
            recdat(rind).Zp{j} = @(phi) interp1(real(recdat(rind).ff{j}), recdat(rind).zzp{j}, phi);
        end
    end
    
    % Flip the path to repeat.
    path = -path;
    
end

%% Plot
figure(1); figshift

% plotme(meshdat.UHP(1).f, real(meshdat.UHP(1).z), 'interp', 10);
hold all
% plotme(meshdat.LHP.f, real(meshdat.LHP.z), 'interp', 10);

for j = 1:length(recdat)
    h(j) = plot3(real(recdat(j).f), imag(recdat(j).f), real(recdat(j).z), 'LineWidth', 4);
    lcell{j} = ['j = ', num2str(j)];
end
legend(h, lcell)
xlabel('phi'); ylabel('psi'); zlabel('real(z)');

%%
figure(2); figshift

% plotme(meshdat.UHP(1).f, imag(meshdat.UHP(1).z), 'interp', 10); 
hold all
% plotme(meshdat.LHP.f, imag(meshdat.LHP.z), 'interp', 10)

for j = 1:length(recdat)
    h(j) = plot3(real(recdat(j).f), imag(recdat(j).f), imag(recdat(j).z), 'LineWidth', 3);
    lcell{j} = ['j = ', num2str(j)];
end

xlabel('phi'); ylabel('psi'); zlabel('imag(z)');
legend(h, lcell)

%%

figure(3); figshift

%[~,zpLHP] = zfree(meshdat.LHP.f, solseries);
%plotme(meshdat.UHP(1).f, real(meshdat.UHP(1).zp), 'interp', 10); 
hold all
%plotme(meshdat.LHP.f, real(zpLHP), 'interp', 10);

for j = 1:length(recdat)
    h(j) = plot3(real(recdat(j).f), imag(recdat(j).f), real(recdat(j).zp), 'LineWidth', 3);
    lcell{j} = ['j = ', num2str(j)];
end

xlabel('phi'); ylabel('psi'); zlabel('real(zp)');
legend(h, lcell)

%%
figure(4); figshift

%[~,zpLHP] = zfree(meshdat.LHP.f, solseries);
%plotme(meshdat.UHP(1).f, imag(meshdat.UHP(1).zp), 'interp', 10);
hold all
%plotme(meshdat.LHP.f, imag(zpLHP), 'interp', 10);

for j = 1:length(recdat)
    h(j) = plot3(real(recdat(j).f), imag(recdat(j).f), imag(recdat(j).zp), 'LineWidth', 3);
    lcell{j} = ['j = ', num2str(j)];
end

xlabel('phi'); ylabel('psi'); zlabel('real(zp)');
legend(h, lcell)

%% Save the rectangular data

name = ['../../stokeswavedata/recdat_steep_rec', num2str(recnum), '_ep', num2str(ep), '.mat'];
save(name, 'recdat');
