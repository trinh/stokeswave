function F = series_bernoulli_finite(a, ep, N, kap, grav_coeff, r0)
% SERIES_BERNOULLI_FINITE defines the residual for Vanden-Broeck's scheme
%
%   This calculates the residual for Bernoulli's equation for the series
%   method as described by Vanden-Broeck. Note that the domain is on 
%   0 <= phi <= 0.5 and assumes symmetry about phi = 0 and periodicity in
%   the band [-0.5, 0.5]
%
%   INPUT
%           a   = N+1 vector of unknowns [series coefficients, mu, beta]
%           ep  = steepness parameter
%           N   = a is of length (N+1)
%           kap = surface tension parameter > 0
%           grav_coeff = allows gravity to be set off
%           r0  = depth as measured in the t-plane
%
%   OUTPUT 
%           F   = N+1 vector of residual values (N Bernoulli + steepness)
%
%   2015-06: Original coding
%   2016-05: Finite-depth version
%   2016-07: Cleaned up

    mu = a(N);
    B = a(N+1);
    
    % Set N+1 zero vector for F
    F = 0*a;    
        
    for j=1:N
        % Set points in phi from 0 to 0.5
        phi = 0.5*((j-1)/(N-1));
        
        % Series expansion (6.13) in JMVB book
        dzdf = 1;
        for k=1:N-1
            t = exp(-2i*pi*phi);
            dzdf = dzdf + a(k)*(t^k + r0^(2*k)*t^(-k));
        end
        
        %The integral term (gravity)
        y = 0;
        for k=1:N-1
            y = y + 1/(2*k*pi)*(cos(2*k*pi*phi)-1)*a(k).*(1 - r0^(2*k));
        end
        
        % ------------------------------
        % Surface tension implementation 
        % ------------------------------
        % for i=1:N-1
        %   d1(k) = -2*pi*a(k)*cos(2*pi*k*phi); %beta_phi
        %   d2(k) = -2*pi*a(k)*sin(2*pi*k*phi); %delta_phi
        % end
        
        % delta = real(dzdf);
        % beta = imag(dzdf);
        % absval_f = (abs(dzdf))^2;
        
        % tension = ((-kap)/(2*pi*mu))*(delta*(sum(d1)) - beta*(sum(d2)))/((absval_f)^(3/2));
        
        tension = 0;
        grav = 2*pi/mu*y;
        energy = 1/(2*abs(dzdf)^2);
        
        %set the jth equation (Bernoulli condition on free surface)
        F(j) = grav_coeff*grav + energy + tension - B;
    end
    
    % For the last equation, use the steepness condition. Vanden-Broeck 
    % (1986, Phys Fluids, 29) suggests using 
    %
    %       ep^2 = 1 - abs(W(0))^2*abs(W(0.5))^2
    %
    % although this is for infinite depth
    
    delta_beta0 = 1 + sum(a(1:N-1).*(1 + r0.^(2*(1:N-1))));
    W0 = 1/(delta_beta0);
    delta_beta_end = 1 + sum(a(1:N-1).* ...
                            (exp(-pi*1i*(1:N-1)) + ...
                             r0.^(2*(1:N-1)).*exp(pi*1i*(1:N-1))));
    Wend = 1/delta_beta_end;
    
    % F(N+1) = 1 - 1/delta0^2 - ep;
    F(N+1) = 1 - abs(W0)^2*abs(Wend)^2 - ep^2;
    
   
end

