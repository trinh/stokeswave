classdef DeepWave
    
    N = 100; % N-1 is no. of terms in power series
% ep = 0.7746; % JMVB 1986 (check mu ~ 1.22)
ep = 0.1;
kap = 0;
grav_coeff = 1;

% Initial guess - last two are mu and B (in that order!)
x0 = zeros(1,N-1);
x0(N) = 1;
x0(N+1) = 0.1;

    
    
    properties (Constant=true)
        kap = 0;            % Surface tension
        grav_coeff = 1;     % Gravity on/off
        mytol = 1e-10;      % Solver tolerance
    end
    properties
        N = 20;             % N-1 is no. of terms in power series
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
            obj.x0 = 0.1*ones(1,obj.N-1);
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
            % ZFREE constructs the series solution using coefficients from solseries
            
            an = obj.an;
            r0 = obj.r0;
            N = obj.N;
            
            z = f;
            zp = 1;
            for k = 1:(N - 1)
                t = exp(-2i*pi*f);
                z = z + an(k)/(2i*pi*k)*(-t.^k + r0^(2*k)*t.^(-k));
                zp = zp + an(k)*(t.^k + r0^(2*k)*t.^(-k));%CHECK IF THIS IS CORRECT
            end
            %     z = z - 1i*solseries.ycrest;
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

