function pathdat = qt_bdint_path_solve(path, wave)
    % Solve ODE along a path
    % s(1) ---- s(2) ---- ... --- s(n-1) --- s(n)
    %
    %   REQUIRED INPUTS
    %       s = path specified as complex vector
    %       wave = wave structure for free surface
    
    % ODE options
   
    mytol = 1e-9;
    options = odeset('RelTol', mytol, 'AbsTol', mytol);
    
    r0 = wave.r0;
    Q = -log(r0)/(2*pi);
    
    L = path.L;
    UD = path.UD;
    
    s = path.f;
    
    if imag(s(1)) > 0
        a = 1;
    else
        a = -1;
    end
    
    pathdat = [];
    sheetnumber = 1;
    
    pathdat.sheet(sheetnumber).f = [];
    pathdat.sheet(sheetnumber).z = [];
    pathdat.sheet(sheetnumber).omega = [];
    pathdat.sheet(sheetnumber).omegap = [];
    pathdat.sheet(sheetnumber).tau = [];
    pathdat.sheet(sheetnumber).taup = [];
    pathdat.sheet(sheetnumber).theta = [];
    
    % For the singularity integral
    integralval = 0;
    
    % Get the last tau and z values
    % Note that dz/df = exp(1i*theta)/q
    [last_z_value, last_dzdf_value] = wave.getZValues(real(s(1)));
    last_q_value = 1/abs(last_dzdf_value);
    last_tau_value = log(last_q_value);
    keyboard
    for j = 1:length(s)-1
        
        z0 = s(j);
        z1 = s(j+1);
        
        if (isnan(z1) == 1)
            % If the next endpoint is NaN
            sheetnumber = sheetnumber + 1; % Increment the sheet
            
            % Supplement initial condition by value of tau(z0) from fluid sheet
            [~, dz0df] = wave.getZValues(real(z0));
            q0 = 1/abs(dz0df);
            last_tau_value = [last_tau_value(:); log(q0)];
              
            % Create new space
            pathdat.sheet(sheetnumber).f = [];
            pathdat.sheet(sheetnumber).tau = [];
            pathdat.sheet(sheetnumber).theta = [];
            pathdat.sheet(sheetnumber).omega = [];
            pathdat.sheet(sheetnumber).omegap = [];
            pathdat.sheet(sheetnumber).z = [];
            continue;
        elseif (isnan(z0) == 1)
            % Do nothing on this step
            continue;
        end
        
        % Compute straight lines.
        
        gam=@(s) z0 + s*(z1-z0);
         
        dgam = (z1 - z0);
        
        % Now solve the ODE with initial conditions
        %fwd = @(s,U)qt_bdint_funccomplex(gam(s), U, wave, sheetnumber, a)*dgam;
        fwd = @(s,U)gen_finite_ODE(gam(s), U, wave, sheetnumber, a, L, UD)*dgam;
        
        %[ss, uu] = ode45(fwd, [0, 1], last_tau_value, options);
        
        tic
        [ss, uu] = ode113(fwd, [0, 1], last_tau_value, options);        
        toc
        % Data for varying ss should be row vectors
        ss = ss.';
        uu = uu.';
        last_tau_value = uu(:, end);
        
        
        % Re-run to get theta
        [taup, tet] = qt_bdint_funccomplex(gam(ss), uu, wave, sheetnumber, a);
        taup = reshape(taup, [sheetnumber, length(ss)]);
        tet = reshape(tet, [sheetnumber, length(ss)]);
        
        pathdat.sheet(sheetnumber).f = [pathdat.sheet(sheetnumber).f, gam(ss)];
        pathdat.sheet(sheetnumber).tau = [pathdat.sheet(sheetnumber).tau, uu];
        pathdat.sheet(sheetnumber).taup = [pathdat.sheet(sheetnumber).taup, taup];
        pathdat.sheet(sheetnumber).theta = [pathdat.sheet(sheetnumber).theta, tet];
        
        % Calculate omega and its derivative
        omega = uu - 1i*tet;
        keyboard
        pathdat.sheet(sheetnumber).omega = [pathdat.sheet(sheetnumber).omega, omega];
        
        % Derivative of omega (in f)
        omegap = 0*omega;
        for k = 1:sheetnumber
            omegap(k,:) = 1/dgam*central_diff(omega(k,:), ss);
        end
        pathdat.sheet(sheetnumber).omegap = [pathdat.sheet(sheetnumber).omegap, omegap];
        
        % Get the z values for the highest sheet
        pathdat.sheet(sheetnumber).z = [pathdat.sheet(sheetnumber).z, ...
            last_z_value + cumtrapz(ss, dgam*exp(-omega(1,:)))];
        last_z_value = pathdat.sheet(sheetnumber).z(end);                
    end
    
    disp(['Requested integral value = ', num2str(integralval, 16)])
    
    % For easy access, compile the top (final) quantities
    pathdat.f = [];
    pathdat.z = [];
    pathdat.tau = [];
    pathdat.theta = [];
    pathdat.omega = [];
    pathdat.omegap = [];
    
    for j = 1:length(pathdat.sheet)
        pathdat.f = [pathdat.f, pathdat.sheet(j).f];
        pathdat.z = [pathdat.z, pathdat.sheet(j).z];
        pathdat.tau = [pathdat.tau, pathdat.sheet(j).tau(1,:)];
        pathdat.theta = [pathdat.theta, pathdat.sheet(j).theta(1,:)];
        pathdat.omega = [pathdat.omega, pathdat.sheet(j).omega(1,:)];
        pathdat.omegap = [pathdat.omegap, pathdat.sheet(j).omegap(1,:)];
        
        
    end
    
%     for j = 2:length(pathdat.sheet)
%         pathdat.tau2 = [pathdat.tau, pathdat.sheet(j).tau(2,:)];
%         pathdat.theta2 = [pathdat.theta, pathdat.sheet(j).theta(2,:)];
%     end
%     
%     for j=3:length(pathdat.sheet)
%         pathdat.tau3 = [pathdat.tau, pathdat.sheet(j).tau(3,:)];
%         pathdat.theta3 = [pathdat.theta, pathdat.sheet(j).theta(3,:)];
%     end