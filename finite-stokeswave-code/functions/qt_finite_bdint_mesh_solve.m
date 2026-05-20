function pathdat = qt_finite_bdint_mesh_solve(s, wave, varargin)
    % Solve ODE along a path
    % s(1) ---- s(2) ---- ... --- s(n-1) --- s(n)
    %
    % This is a modification of the path_solve code. It first solves along a
    % path to a predetermined location (marked by 999) then it meshes either
    % the UHP or LHP (depending on where the location is set)
    %
    % 2026-05-20: Modification by Phil to get the finite depth version
    % working; changes include making sure that we are leveraing the use of
    % the VBStokesWave class instead of individual functions 
    %
    % Input: 
    %   s = a path specification via get_path()
    %   wave = a VBStokesWave object
    
    % ODE options
    mytol = 1e-9;
    options = odeset('RelTol', mytol, 'AbsTol', mytol);
    
    % Parse the inputs
    myparse = inputParser;
    myparse.addRequired('path', @(x) true);
    myparse.addRequired('wave', @(x) true);    
    myparse.addParameter('PsiMax', 1, @(x) x > 0)
    myparse.addParameter('PsiMin', -1, @(x) x < 0)
    myparse.addParameter('NPsi', 30, @(x) x > 0)
    myparse.addParameter('NPhi', 30, @(x) x > 0)
    myparse.parse(path, wave, varargin{:})
       
    if imag(s(1)) > 0
        a = 1;
    else
        a = -1;
    end
    
    pathdat = [];
    sheetnumber = 1;
    
    pathdat.sheet(sheetnumber).f = [];
    pathdat.sheet(sheetnumber).omega = [];
    pathdat.sheet(sheetnumber).omegap = [];
    pathdat.sheet(sheetnumber).tau = [];
    pathdat.sheet(sheetnumber).taup = [];
    pathdat.sheet(sheetnumber).theta = [];
    
    % Create the meshing
    psimax = myparse.Results.PsiMax;
    psimin = myparse.Results.PsiMin;
    N_psi = myparse.Results.NPsi;
    N_phi = myparse.Results.NPhi;
    
    % No initial value was specified, use the wave object
    [last_q_value, ~] = wave.getQTValues(real(s(1)));
    [last_z_value, ~] = wave.getZValues(real(s(1)));
    last_tau_value = log(last_q_value);
        % last_tau_value = log(get_q_series(real(s(1)), solseries));
    % [last_z_value, ~] = zfree(real(s(1)), solseries);
    
    for j = 1:length(s)-1
        
        z0 = s(j);
        z1 = s(j+1);
        
        if (isnan(z1) == 1)
            % If the next endpoint is NaN
            sheetnumber = sheetnumber + 1; % Increment the sheet
            
            % Supplement initial condition by value of tau(z0) from fluid sheet
            [q_tmp, ~] = wave.getQTValues(real(z0));
            last_tau_value = [last_tau_value(:); log(q_tmp)];
            
            % Create new space
            pathdat.sheet(sheetnumber).f = [];
            pathdat.sheet(sheetnumber).tau = [];
            pathdat.sheet(sheetnumber).theta = [];
            continue;
        elseif (isnan(z0) == 1)
            % Do nothing on this step
            continue;
        elseif z1 == 999
            % This marks when we want to begin the meshing
            phi0 = abs(real(z0));
            psi0 = imag(z0);
            phi = linspace(-phi0, phi0, N_phi);
            if psi0 < 0
                psi = linspace(psi0, psimin, N_psi);
            else
                psi = linspace(psi0, psimax, N_psi);
            end
            [PHI, PSI] = meshgrid(phi, psi);
            FF = PHI + 1i*PSI;
            taumat = 0*FF;
            thetamat = 0*FF;
            zmat = 0*FF;
            
            % First solve across
            fwd = @(s,U)qt_finite_bdint_funccomplex(s + 1i*psi0, U,  wave, sheetnumber, a);
            %         [ss, taurow] = ode45(fwd, phi, last_tau_value, options);
            [ss, taurow] = ode113(fwd, phi, last_tau_value, options);
            thetarow = 0*taurow;
            for kk = 1:length(ss)
                % Have to convert the input back to a column
                tmp = taurow(kk,:);
                tmp = tmp(:);
                [~, thetaval] = qt_finite_bdint_funccomplex(ss(kk) + 1i*psi0, tmp, wave, sheetnumber, a);
                % Only record the top value
                thetarow(kk,:) = thetaval;
            end
            
            omegarow = taurow(:,1) - 1i*thetarow(:,1);
            zrow = last_z_value + cumtrapz(ss, exp(-omegarow));
            
            % Now solve in vertical
            for k = 1:length(phi)
                fprintf('%d ', k);
                disp(['phi(k) = ', num2str(phi(k))]);
                
                fwd = @(s,U) 1i*qt_finite_bdint_funccomplex(phi(k) + 1i*s, U, wave, sheetnumber, a);
                %             [ss, uu] = ode45(fwd, psi, taurow(k,:));
                [ss, uu] = ode113(fwd, psi, taurow(k,:), options);
                % Re-run to get theta
                tt = 0*uu;
                for kk = 1:length(ss)
                    % Have to convert input back to a column
                    tmp = uu(kk,:);
                    tmp = tmp(:);
                    [~, thetaval] = qt_finite_bdint_funccomplex(phi(k) + 1i*ss(kk), tmp, wave, sheetnumber, a);
                    tt(kk,:) = thetaval;
                end
                
                taumat(:,k) = uu(:,1);
                thetamat(:,k) = tt(:,1);
                
                % We integrate to get the value of z
                zmat(:,k) = zrow(k) + cumtrapz(ss, 1i*exp(-(taumat(:,k) - 1i*thetamat(:,k))));                
            end
            
            
            
            % Flip if necessary
            if psi0 < 0
                taumat = flipud(taumat);
                thetamat = flipud(thetamat);
                zmat = flipud(zmat);
                FF = flipud(FF);
            end
            
            pathdat.taumat = taumat;
            pathdat.thetamat = thetamat;
            pathdat.FF = FF;
            pathdat.zmat = zmat;
            break;
            
        end
        
        % Compute straight line
        gam = @(s) z0 + s*(z1 - z0);
        dgam = (z1 - z0);
        
        % Now solve the ODE with initial conditions
        fwd = @(s,U)qt_finite_bdint_funccomplex(gam(s), U, wave, sheetnumber, a)*dgam;
        %     [ss, uu] = ode45(fwd, [0, 1], last_tau_value, options);
        [ss, uu] = ode113(fwd, [0, 1], last_tau_value, options);
        
        % Data for varying ss should be row vectors
        ss = ss.';
        uu = uu.';
        last_tau_value = uu(:, end);
        
        % NOTE: the sheetnumber starts at 1 for the desired solution with each
        % addition + 1 for the residue
        
        % Re-run to get theta
        [taup, tet] = qt_finite_bdint_funccomplex(gam(ss), uu, wave, sheetnumber, a);
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
        
        % Get the value of z at the end
        % Note here we use the sheetnumber = 1
        last_z_value = last_z_value + trapz(ss, dgam*exp(-omega(1,:)));
        
        
    end
    % keyboard
    
    % For easy access, compile the top (final) quantities
    pathdat.f = [];
    pathdat.tau = [];
    pathdat.theta = [];
    pathdat.omega = [];
    pathdat.omegap = [];
    
    for j = 1:length(pathdat.sheet)
        %     j
        pathdat.f = [pathdat.f, pathdat.sheet(j).f];
        pathdat.tau = [pathdat.tau, pathdat.sheet(j).tau(1,:)];
        pathdat.theta = [pathdat.theta, pathdat.sheet(j).theta(1,:)];
        pathdat.omega = [pathdat.omega, pathdat.sheet(j).omega(1,:)];
        pathdat.omegap = [pathdat.omegap, pathdat.sheet(j).omegap(1,:)];
    end
