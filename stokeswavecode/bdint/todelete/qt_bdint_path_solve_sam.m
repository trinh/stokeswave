function [pathdat, singloc] = qt_bdint_path_solve_sam(s, solseries, sing_start)
% Solve ODE along a path 
% s(1) ---- s(2) ---- ... --- s(n-1) --- s(n)

if imag(s(1)) > 0
    a = 1;
else
    a = -1;
end

pathdat = [];
sheetnumber = 1;

pathdat.sheet(sheetnumber).f = [];
pathdat.sheet(sheetnumber).tau = [];
pathdat.sheet(sheetnumber).theta = [];

% For the singularity integral
integralval = 0;

last_tau_value = log(get_q_series(real(s(1)), solseries));

for j = 1:length(s)-1
    
    z0 = s(j);
    z1 = s(j+1);
                  
    if (isnan(z1) == 1)
        % If the next endpoint is NaN
        sheetnumber = sheetnumber + 1; % Increment the sheet
        
        % Supplement initial condition by value of tau(z0) from fluid sheet
        last_tau_value = [last_tau_value(:); log(get_q_series(real(z0), solseries))];
        
        % Create new space
        pathdat.sheet(sheetnumber).f = [];
        pathdat.sheet(sheetnumber).tau = [];
        pathdat.sheet(sheetnumber).theta = [];
        continue;
    elseif (isnan(z0) == 1)
        % Do nothing on this step
        continue;
    end
    
    % Compute straight line
    gam = @(s) z0 + s*(z1 - z0);
    dgam = (z1 - z0);
    
    % Now solve the ODE with initial conditions
    fwd = @(s,U)qt_bdint_funccomplex(gam(s), U, solseries, sheetnumber, a)*dgam;   
    [ss, uu] = ode45(fwd, [0, 1], last_tau_value);

    % Data for varying ss should be row vectors
    ss = ss.';
    uu = uu.';
    last_tau_value = uu(:, end);
    
    
    % Re-run to get theta
    [taup, tet] = qt_bdint_funccomplex(gam(ss), uu, solseries, sheetnumber, a);
    taup = reshape(taup, [sheetnumber, length(ss)]);
    tet = reshape(tet, [sheetnumber, length(ss)]);
    
    pathdat.sheet(sheetnumber).f = [pathdat.sheet(sheetnumber).f, gam(ss)];
    pathdat.sheet(sheetnumber).tau = [pathdat.sheet(sheetnumber).tau, uu];  
    pathdat.sheet(sheetnumber).theta = [pathdat.sheet(sheetnumber).theta, tet];   
    
    
    if j >= sing_start
    % Calculate the integral for the singlarity location
    zzp = 1./(exp(uu(1,:) - 1i.*tet(1,:)));
    zzpp = central_diff(zzp, ss);
    
    integrand = gam(ss).*zzpp./zzp;
    integralval = integralval + trapz(ss,integrand);
    
    hold on
    plot3(real(gam(ss)), imag(gam(ss)), imag(zzpp));
    end
    
   
end

singloc = (1/(2*pi*1i))*integralval;

% For easy access, compile the top (final) quantities
pathdat.f = [];
pathdat.tau = [];
pathdat.theta = [];
for j = 1:length(pathdat.sheet)
    pathdat.f = [pathdat.f, pathdat.sheet(j).f];
    pathdat.tau = [pathdat.tau, pathdat.sheet(j).tau(1,:)];
    pathdat.theta = [pathdat.theta, pathdat.sheet(j).theta(1,:)];
end
