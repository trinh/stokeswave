clear
close all

% =====================================
%         Load in the solution
F = load('data/series_N200_ep0.9.mat');
solseries = F.solseries;

F = load('data/meshdat_M100_ep0.9.mat');
meshdat = F.meshdat;

mu = solseries.mu;
B = solseries.B;
% =====================================


mytol = 1e-7;
options = odeset('RelTol', mytol, 'AbsTol', mytol);

% This is the main data structure
recdat = [];

phi0 = 0.4;
psi0 = 1;

path = [phi0, phi0 - 1i*psi0, -phi0 - 1i*psi0, -phi0, phi0];

M = 50; % 1/4 number of points in the in-fluid mesh


for rind = 1:15  
    recdat(rind).f = [];
    recdat(rind).z = [];
    recdat(rind).zp = [];
    for j = 1:length(path)-1
        
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
            if j == 3
%                 keyboard
            end
        else
            % Left/right
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
            
            if j == 1
                z0 = recdat(rind-1).zz{3}(end);
            else
                z0 = zz(end);
            end
            
            fwd = @(s, z) dgam*zbern_flipode(s, z, solseries, Zminus, Zpminus);
                    
            
            [ss, zz] = ode45(fwd, [0, 1], z0, options);
            
            zzp = zbern_flipode(ss, zz, solseries, Zminus, Zpminus);
            zzp = zzp(:);
            ss = ss(:);
            zz = zz(:);
            
            recdat(rind).ff{j} = gam(ss);
            recdat(rind).zz{j} = zz;
            recdat(rind).zzp{j} = zzp;
        end
        recdat(rind).f = [recdat(rind).f; recdat(rind).ff{j}];
        recdat(rind).z = [recdat(rind).z; recdat(rind).zz{j}];
        recdat(rind).zp = [recdat(rind).zp; recdat(rind).zzp{j}];
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
    
    path = -path;
    
end

%%
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

xlabel('phi'); ylabel('psi'); imag('real(z)');
legend(h, lcell)

%%

figure(3); figshift

[~,zpLHP] = zfree(meshdat.LHP.f, solseries);
plotme(meshdat.UHP(1).f, real(meshdat.UHP(1).zp), 'interp', 10); 
hold all
plotme(meshdat.LHP.f, real(zpLHP), 'interp', 10);

for j = 1:length(recdat)
    h(j) = plot3(real(recdat(j).f), imag(recdat(j).f), real(recdat(j).zp), 'LineWidth', 3);
    lcell{j} = ['j = ', num2str(j)];
end

xlabel('phi'); ylabel('psi'); zlabel('real(zp)');
legend(h, lcell)

%%
figure(4); figshift

[~,zpLHP] = zfree(meshdat.LHP.f, solseries);
plotme(meshdat.UHP(1).f, imag(meshdat.UHP(1).zp), 'interp', 10) 
hold all
plotme(meshdat.LHP.f, imag(zpLHP), 'interp', 10);

for j = 1:length(recdat)
    h(j) = plot3(real(recdat(j).f), imag(recdat(j).f), imag(recdat(j).zp), 'LineWidth', 3);
    lcell{j} = ['j = ', num2str(j)];
end

xlabel('phi'); ylabel('psi'); zlabel('real(zp)');
legend(h, lcell)
