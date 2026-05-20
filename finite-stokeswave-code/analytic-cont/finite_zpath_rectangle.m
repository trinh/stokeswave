clear
  clear classes
% close all

% =====================================
%         Load in the solution
r0 = 0.3;
ep = sqrt(0.1);
wave = VBStokesWave(r0, ep);
% =====================================

%psifloor = 1/(2*pi)*log(r0);

mytol = 1e-9;
options = odeset('RelTol', mytol, 'AbsTol', mytol);



% This is the main data structure
recdat = [];

phi0 = 0.4;
psi0 = 0.1;

path = [phi0, phi0 + 1i*psi0, -phi0 + 1i*psi0, -phi0,phi0];

M = 50; % 1/4 number of points in the in-fluid mesh


for rind = 1:2
    
    recdat(rind).f = [];
    recdat(rind).z = [];
    recdat(rind).zp = [];
    
    if rind == 1 %only plot psi = 0 path on physical domain
        lpath = length(path)-1;
    else lpath = length(path)-2;
    end
    
    for j = 1:lpath
        
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
            [recdat(rind).zz{j}, recdat(rind).zzp{j}] = wave.getZValues(gam(ss));                        
            
        else
            % Left/right
            if j == 1 || j == 3
                % Set up a function from (rind - 1) data and evaluate
                % previous interpolation given a value for s
                Zminus = @(s) recdat(rind-1).Z{j}(-imag(gam(s)));
                Zpminus = @(s) recdat(rind-1).Zp{j}(-imag(gam(s)));                                
            elseif j == 2 || j == 4
                Zminus = @(s) recdat(rind-1).Z{j}(-real(gam(s)));
                Zpminus = @(s) recdat(rind-1).Zp{j}(-real(gam(s)));
            end
            
            % Initial condition for ODE
            if j == 1
                % If starting a new rectangle, use the end of the 3rd path
                % (LHS of old rectangle)
                z0 = recdat(rind-1).zz{3}(end);
            else
                % If continuing rectangle, use the previous z-segment-value
                z0 = zz(end);
            end
            
            % Create the ODE function file
            % dz/df = zbern_flipode
            % dz/ds = (dgamma/ds)*dzdf
            fwd = @(s, z) dgam*finite_zbern_flipode(s, z, wave, Zminus, Zpminus);
                    
            % Solve the ODE along s in [0, 1]
            [ss, zz] = ode45(fwd, [0, 1], z0, options);
            
            % Retrieve dz/df values
            zzp = finite_zbern_flipode(ss, zz, wave, Zminus, Zpminus);
            zzp = zzp(:);
            ss = ss(:);
            zz = zz(:);
            
            recdat(rind).ff{j} = gam(ss);
            recdat(rind).zz{j} = zz;
            recdat(rind).zzp{j} = zzp;
        end
        
        % Join the new data with the old (for plotting)
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
    
    % Negate the path in preparation for next
    path = -path;    
end


name = ['../../finite-stokeswave-data/recdata_ep', num2str(ep), '_r0', num2str(r0), '.mat'];
save(name);


% Plot real(z)

figure(2); figshift

%plotme(meshdat.UHP(1).f, real(meshdat.UHP(1).z), 'interp', 10);
hold all
%plotme(meshdat.LHP.f, real(meshdat.LHP.z), 'interp', 10);

for j = 1:length(recdat)
    h(j) = plot3(real(recdat(j).f), imag(recdat(j).f), real(recdat(j).z), 'LineWidth', 4);
    lcell{j} = ['j = ', num2str(j)];
end
legend(h, lcell)
xlabel('phi'); ylabel('psi'); zlabel('real(z)');
view([-40, 30]);


%% Plot imag(z)

figure(1); figshift

% plotme(meshdat.UHP(1).f, imag(meshdat.UHP(1).z), 'interp', 10); 
hold all
% plotme(meshdat.LHP.f, imag(meshdat.LHP.z), 'interp', 10)

for j = 1:length(recdat)
    h(j) = plot3(real(recdat(j).f), imag(recdat(j).f), imag(recdat(j).z), 'LineWidth', 3);    
    lcell{j} = ['j = ', num2str(j)];
end

%line([-phi0 phi0],[psifloor psifloor]);

% mesh(z);
% % Z threshold value. 
% threshold = .5; % please change this as needed .
% % Obtain the limits of the axes
% yp = get(gca,'Ylim');
% xp = get(gca,'Xlim');
% % Use the axes x and Y limits to find the co-ordinates for the patch
% phi1 = [ -phi0 phi0];
% z1 = [-max(imag(recdat(j).z)), max(imag(recdat(j).z))];
% 
% psi1 = psifloor*ones(1,numel(x1))* threshold;  % creates a 1x4 vector representing the Z coordinate value
% 
% p = patch(x1,y1,z1, 'b');
% 
% % Set the Face and edge transparency to 0.2 using the following properties
% set(p,'facealpha',0.2)
% set(p,'edgealpha',0.2)

xlabel('phi'); ylabel('psi'); imag('real(z)');
legend(h, lcell)
view([-40, 30]);

