function pathdat = qt_bdint_path_solve(path, solseries, varargin)
    % Solve ODE along a path
    % s(1) ---- s(2) ---- ... --- s(n-1) --- s(n)
    %
    %   REQUIRED INPUTS
    %       s = path specified as complex vector
    %       solseries = solution structure
    %   OPTIONAL INPUT (NAME/VALUE PAIRS)
    %       'singcoef' with b where b is the singularity location
    %
    %       e.g. qt_bdint_path_solve(s, solseries, 'singcoef', 0.15)
    
    % ODE options
    mytol = 1e-9;
    options = odeset('RelTol', mytol, 'AbsTol', mytol);
    
    % Parse the inputs
    myparse = inputParser;
    myparse.addRequired('pathdat',@(x) true);
    myparse.addRequired('series', @isstruct);
    myparse.addParameter('singcoef', NaN, @(x) isnumeric);
    myparse.parse(path, solseries, varargin{:})
    
    s = path.f;
    singstart = path.singstart;
    
    % If we had requested to track singularity location or coefficient
    if isnan(myparse.Results.singcoef) == 0
        singloc = myparse.Results.singcoef;
        % Check that singstart is defined
        if singstart < Inf
            tracktype = 'singcoef';
        else
            warning('something is wrong in singcoef');
            keyboard
        end
    elseif singstart < Inf
        tracktype = 'singloc';
    else
        tracktype = 'none';
    end
    
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
    
    last_tau_value = log(get_q_series(real(s(1)), solseries));
    [last_z_value, ~] = zfree(real(s(1)), solseries);
    
    
    for j = 1:length(s)-1
        
        z0 = s(j);
        z1 = s(j+1);
        
        if (isnan(z1) == 1)
            % If the next endpoint is NaN
            sheetnumber = sheetnumber + 1; % Increment the sheet
            
            % Supplement initial condition by value of tau(z0) from fluid sheet
            % keyboard
            last_tau_value = [last_tau_value(:); log(get_q_series(real(z0), solseries))];
            
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
        
        % Compute straight line
        gam = @(s) z0 + s*(z1 - z0);
        dgam = (z1 - z0);
        
        % Now solve the ODE with initial conditions
        fwd = @(s,U)qt_bdint_funccomplex(gam(s), U, solseries, sheetnumber, a)*dgam;
        %         [ss, uu] = ode45(fwd, [0, 1], last_tau_value, options);
        tic
        [ss, uu] = ode113(fwd, [0, 1], last_tau_value, options);
        toc
        
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
        pathdat.sheet(sheetnumber).taup = [pathdat.sheet(sheetnumber).taup, taup];
        pathdat.sheet(sheetnumber).theta = [pathdat.sheet(sheetnumber).theta, tet];
        
        % Calculate omega and its derivative
        omega = uu - 1i*tet;
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
        
        if j >= singstart
            % SAM: you may be able to just use my omegap function here but I've
            % kept yours here unaltered
            % Calculate the integral for the singlarity location
            zzp = 1./(exp(uu(1,:) - 1i.*tet(1,:)));
            zzpp = central_diff(zzp, ss);
            
            switch tracktype
                case 'singloc'
                    integrand = gam(ss).*zzpp./zzp;
                    integralval = integralval + 1/(2*pi*1i)*trapz(ss,integrand);
                case 'singcoef'
                    % integrand = 2*dgam.*zzp.*(gam(ss)-singloc).^(-1/2); % Check this code
                    integrand = dgam*(gam(ss)-singloc).*zzpp.*zzp;
                    integralval = integralval + 1/(2*pi*1i)*trapz(ss,integrand);
            end
            
        end
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
    
    
