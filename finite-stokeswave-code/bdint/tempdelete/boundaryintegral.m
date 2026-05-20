function F = boundaryintegral(X, N, ep, r0)

%X is vector of unknowns -- 2N size for N delta and N beta
%N should be odd - mesh size
%ep is the wave steepness parameter

% Unpack the solution vector
delta = X(1:N);
beta = 0*delta;
beta(1) = 0;
beta(N) = 0;
beta(2:N-1) = X(N+1:2*N-2);    
B = X(2*N-1);
mu = X(2*N);

%mesh points
phi = linspace(-0.5,0.5,N);
h = phi(2) - phi(1);

%midpoints 
phi_m =  1/2*(phi(1:N-1) + phi(2:N));

%interpolated midpoint values
% delta_m = 1/2*(delta(1:N-1) + delta(2:N));
% beta_m = 1/2*(beta(1:N-1) + beta(2:N));
%4 pt interpolation (p.154 JMVB)
% delta_m = 1/16*[(-delta(N-1) + 9*delta(1) + 9*delta(2) - delta(3)) ...
%                 (-delta(1:N-3) + 9*delta(2:N-2) + 9*delta(3:N-1) - delta(4:N)) ...
%                 (-delta(N-2) + 9*delta(N-1) + 9*delta(N) - delta(2))];
% beta_m = 1/16*[(-beta(N-1) + 9*beta(1) + 9*beta(2) - beta(3)) ...
%                 (-beta(1:N-3) + 9*beta(2:N-2) + 9*beta(3:N-1) - beta(4:N)) ...
%                 (-beta(N-2) + 9*beta(N-1) + 9*beta(N) - beta(2))];

% Interpolate delta and beta at midpoints using Matlab spline
delta_m = interp1(phi, delta, phi_m, 'spline');
beta_m = interp1(phi, beta, phi_m, 'spline');

% Error vectors
F = zeros(2*N, 1);

%evaluate the boundary integral terms at each point (first N-1 eqns)
for k=1:N-1
    % Integrand of integral
    y = 2*real((delta + 1i*beta).*exp(-2i*pi*phi)./(exp(-2i*pi*phi) - exp(-2i*pi*phi_m(k)))) ...
        - 2*r0^2*real((delta + 1i*beta)*r0^2.*exp(-2i*pi*phi)./(r0^2*exp(-2i*pi*phi) - exp(-2i*pi*phi_m(k))));
    
    % Significant speedup if we code our own trapezoid but if doubt, verify
    % with canned routine    
    g = 2*y;
    g(1) = g(1)/2;
    g(N) = g(N)/2;
    integral = h/2*sum(g);    
    F(k) = delta_m(k) - integral;
end

%now we evaluate the Bernoulli equation

%first the integral term.
m=(N+1)/2;

% Y = cumtrapz(phi, beta);
% Y = Y - Y(m);
Y = 0*delta;
for j = m:N-1
    Y(j+1) = Y(j) + h/2*(beta(j+1) + beta(j));
end
for j = m:-1:2
    Y(j-1) = Y(j) - h/2*(beta(j-1) + beta(j));
end

% Get midpoints for Y
% Y_m = 1/2*(Y(1:N-1) + Y(2:N));

% 4-point interolation and use periodicity to handle the first and last
% points
Y_m = 1/16*[(-Y(N-1) + 9*Y(1) + 9*Y(2) - Y(3)) ...
                (-Y(1:N-3) + 9*Y(2:N-2) + 9*Y(3:N-1) - Y(4:N)) ...
                (-Y(N-2) + 9*Y(N-1) + 9*Y(N) - Y(2))];

% Bernoulli's equation for N-1 equations
F(N:2*N-2) = 1./(delta_m.^2 + beta_m.^2) + ((2*pi)/mu)*Y_m - B;
% Using the below seems to get rid of the sawtooth, but it gives slight
% asymetry to the wave
% boo = 1./(delta.^2 + beta.^2) + ((2*pi)/mu)*Y - B;
% F(N:2*N-2) = boo(1:N-1);

%equations for B and mu. (fixing wavelength and steepness)
F(2*N-1) = trapz(phi,delta)-1;
%F(2*N) = Y(end) - s;

% Wave amplitude
F(2*N) = 1 - ((1)/(delta(1)^2*delta(m)^2)) - ep^2;

% keyboard
