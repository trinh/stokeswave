function F = boundaryintegral(X, phi, ep)

%X is the evaluation vector
%N should be odd - mesh size
%ep is the wave steepness parameter

N = length(phi);

%Solution vector
delta = X(1:N);
beta = 0*delta;
beta(1) = 0;
beta(N) = 0;
beta(2:N-1) = X(N+1:2*N-2);  

mu = X(2*N-1);
B = X(2*N);

%mesh points
h = phi(2) - phi(1);

%midpoints 
phi_m =  1/2*(phi(1:N-1) + phi(2:N));

% Interpolation at midpoints
%   The scheme we use can affect the sawtooth patterns

% ** Spline midpoints ** : no sawtooth
% delta_m = interp1(phi, delta, phi_m, 'spline');
% beta_m = interp1(phi, beta, phi_m, 'spline');
%
% Periodic 3-point asymmetric interpolation 
D = midpoint_matrix(N);
delta_m = D*delta(:); delta_m = delta_m.';
beta_m = D*beta(:); beta_m = beta_m.';

F = zeros(2*N, 1);

%evaluate the boundary integral terms at each point (first N-1 eqns)
for k=1:N-1
    
    f = beta.*cot(pi*(phi - phi_m(k)));
    % Significant speedup if we code our own trapezoid but if doubt, verify
    % with canned routine
    % integral = trapz(phi,y);
    g = 2*f;
    g(1) = g(1)/2;
    g(N) = g(N)/2;
    integral = h/2*sum(g);
    F(k) = delta_m(k) + integral - 1;
end

%now we evaluate the Bernoulli equation

%first the integral term.
m=(N+1)/2;

Y = 0*delta;
for j = m:N-1
    Y(j+1) = Y(j) + h/2*(beta(j+1) + beta(j));
end
for j = m:-1:2
    Y(j-1) = Y(j) - h/2*(beta(j-1) + beta(j));
end
Y_m = D*Y(:); Y_m = Y_m.';

% Bernoulli's
F(N:2*N-2) = 1/2*1./(delta_m.^2 + beta_m.^2) + ((2*pi)/mu)*Y_m - B;

%equations for B and mu. (fixing wavelength and steepness)
F(2*N-1) = trapz(phi,delta)-1;
%F(2*N) = Y(end) - s;

W0 = 1/(delta(m) + 1i*beta(m));
Wend = 1/(delta(end) + 1i*beta(end));
F(2*N) = 1 - abs(W0)^2*abs(Wend)^2 - ep^2;


function D = midpoint_matrix(N)
    e = ones(N,1);
    
    % ** 3-point asymmetric interpolation 
    % from Fornberg "Generation of finite difference formulas" p. 705
    % This creates a midpoint operator matrix D
    D = spdiags(1/8*[3*e 6*e -e], [0 1 2], N-1, N);
    D(N-1, 2) = -1/8; % periodic; only works if 2 = N-1
    %
    % ** 2-point centred midpoint
    % D = spdiags([e/2  e/2], [0 1], N-1, N);

