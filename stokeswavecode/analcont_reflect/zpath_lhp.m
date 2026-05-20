% Code to solve on an arbitrary path in the UHP (the sheet starting above
% the fluid)

clear
close all


% =====================================
%         Load in the solution
F = load('data/series_N100_ep0.05.mat');
solseries = F.solseries;

mu = solseries.mu;
B = solseries.B;
ep = solseries.ep;
% =====================================

% Path starts on the axis of the first unphysical "fluid" sheet
path = [];

% First we solve on the same reverse path in the UHP of the fluid so we can
% use this data later.

path_minus = -path;
reverse_path_dat = [];

for j=1:length(path)-1       
    % Path segment end points
    f0 = path_minus(j);
    f1 = path_minus(j+1);

    % Compute straight line
    gam = @(s) f0 + s*(f1 - f0);
    dgam = (f1 - f0);
    
    
    [reverse_dat(j).ff, reverse_dat(j).gg, reverse_dat(j).ggp, reverse_dat(j).zzp] = gbern_contour(s,solseries,g0);
end

%{
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
        else
            % Left/right (getting z(-f) and z'(-f) along the previous path segment)
            if j == 1 || j == 3
                Zminus = @(s) recdat(rind-1).Z{j}(-imag(gam(s)));
                Zpminus = @(s) recdat(rind-1).Zp{j}(-imag(gam(s)));
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
            
            zzp = zbern_flipode(ss, zz, solseries, Zminus, Zpminus);
            zzp = zzp(:);
            ss = ss(:);
            zz = zz(:);
            
            % Save the rectangle data
            recdat(rind).ff{j} = gam(ss);
            recdat(rind).zz{j} = zz;
            recdat(rind).zzp{j} = zzp;
        end
        
        % Combine the rectangle data along path segments
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
    
    % Flip the path to repeat.
    path = -path;
%}
    


