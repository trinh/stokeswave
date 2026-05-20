classdef VBStokesWave
    % VBSTOKESWAVE creates the class file for generating a nonlinear Stokes
    % wave in finite depth 
    %
    % Get the coefficients of z'(f) via series truncation and collocation
    %    Finite depth version
    % -----------------------------
    %
    %   2015-06: Original infinite-depth version
    %   2016-05: Finite-depth version & Class file
    %   2016-07: Clean-up 
    
    properties (Constant=true)
        kap = 0;            % Surface tension
        grav_coeff = 1;     % Gravity on/off
        mytol = 1e-13;      % Solver tolerance
    end
    properties
        N = 50;             % N-1 is no. of terms in power series
        % ep = 0.7746;      % JMVB 1986 (check mu ~ 1.22)
        ep;                 % Steepness parameter
        r0;                 % t-plane value for bottom        
        
        x0;                 % Initial guess
        
        options;            % Solver options
        
        % Solver output
        an;
        mu;
        B;
        ycrest;
        H;
    end
    
    methods
        function obj = VBStokesWave(r0, ep)
            % Initial guess - last two are mu and B (in that order!)
            obj.x0 = 0*ones(1,obj.N-1);
            obj.x0(obj.N) = 1;
            obj.x0(obj.N+1) = 0.1;
            
            obj.r0 = r0;
            obj.ep = ep;
            
            obj.options = optimset('Display', 'iter', 'TolFun', obj.mytol, 'TolX', obj.mytol);
            
            obj = getCoef(obj);
        end
        function obj = getCoef(obj)
            % Solve the ODE
            N = obj.N;
            fwd = @(x)series_bernoulli_finite(x, obj.ep, obj.N, obj.kap, obj.grav_coeff, obj.r0);
            [sol,fval] = fsolve(fwd, obj.x0, obj.options); % Call solver
            
            obj.an = sol(1:N-1);
            obj.mu = sol(N);
            obj.B = sol(N+1);
            obj.ycrest = sum(sol(1:N-1)./(2*pi*(1:N-1)));
            obj.H = imag(getZValues(obj, 0)-getZValues(obj, 0.5));
        end
        function [z, zp] = getZValues(obj, f)
            % getZValues constructs the series solution using coefficients from solseries
            
            an = obj.an;
            r0 = obj.r0;
            N = obj.N;
            
            z = f;
            zp = 1;
            for k = 1:(N - 1)
                z = z + an(k)/(2i*pi*k)* ...
                    ( -exp(-2i*pi*f*k) + r0^(2*k)*exp(2i*pi*f*k) ...
                            + (1 - r0^(2*k)) ...
                    );
                zp = zp + an(k)*(exp(-2i*pi*f*k) + r0^(2*k)*exp(2i*pi*f*k));
            end
            
            % No longer required with the modification to z
            % z = z - 1i*obj.ycrest;
        end
        function [q, theta] = getQTValues(obj, f)
            [~, zp] = obj.getZValues(f);
            q = 1./abs(zp);
            theta = atan(imag(zp)./real(zp));    
        end
        function plotStreamlines(obj)
            % Plot the solution
            
            phi = linspace(-0.5,0.5,40);
            if obj.r0 ~= 0
                psi = linspace(1/(2*pi)*log(obj.r0), 0, 40);
            else
                psi = linspace(-2, 0, 40);
            end
            [PHI, PSI] = meshgrid(phi, psi);
            F = PHI + 1i*PSI;
            
            [Z, ZP] = getZValues(obj, F);
            
            figure(1); clf(1); hold all
            for j = 1:length(psi)
                plot(Z(j,:));
            end
        end
    end

end

