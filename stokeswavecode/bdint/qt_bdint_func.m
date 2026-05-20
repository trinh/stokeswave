function F = qt_bdint_func(UU, phi, ep)
% Boundary integral modification to the (q, theta) plane
% This is used to solve for the free surface using boundary integrals

% X is the evaluation vector
% N should be odd - mesh size
% ep is the wave steepness parameter

N = length(phi);

%Solution vector
q = UU(1:N);
theta = 0*q;
theta(1) = 0;
theta(N) = 0;
theta(2:N-1) = UU(N+1:2*N-2);  

mu = UU(2*N-1);
B = UU(2*N);

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
q_m = D*q(:); q_m = q_m.';
theta_m = D*theta(:); theta_m = theta_m.';

%evaluate the boundary integral terms at each point (first N-1 eqns)
bdintrhs_m = 0*phi_m;
for k=1:N-1
    f = theta.*cot(pi*(phi - phi_m(k)));
    % Significant speedup if we code our own trapezoid but if doubt, verify
    % with canned routine
    % integral = trapz(phi,y);
    g = 2*f;
    g(1) = g(1)/2;
    g(N) = g(N)/2;
    bdintrhs_m(k) = h/2*sum(g);    
end

%now we evaluate the Bernoulli equation

%first the integral term.
m=(N+1)/2;
[~, Y] = qt_getxy(phi, q, theta);
Y_m = D*Y(:); Y_m = Y_m.';


F = zeros(2*N, 1);
% Integral equation
F(1:N-1) = log(q_m) - bdintrhs_m; 

% Bernoulli's
% Need derivative
% dqdphi_m = central_diff(q_m, phi_m);
% F(N:2*N-2) = q_m.^2.*dqdphi_m + 2*pi/mu*sin(theta_m);
F(N:2*N-2) = 1/2*q_m.^2 + 2*pi/mu*Y_m - B;

%equations for B and mu. (fixing wavelength and steepness)
delta = cos(theta)./q;
F(2*N-1) = trapz(phi,delta)-1;

% Steepness condition
F(2*N) = 1 - q(m)^2*q(end)^2 - ep^2;


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

