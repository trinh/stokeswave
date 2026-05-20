function F = qt_finite_bdint_func(UU, phi, ep, wave)
% Boundary integral modification to the (q, theta) plane
% This is used to solve for the free surface using boundary integrals

% X is the evaluation vector
% N should be odd - mesh size
% ep is the wave steepness parameter
%
%          May 2026
%  ---------------------------------
%  To help track down possible issues, we are re-writing a boundary
%  integral checker for the finite-depth case. This was taken from
%  qt_bdint_func.m code that was previously working for the infinite-depth
%  Stokes wave 
%
%  For inspiration, look at qt_bdint_funccomplex.m in the finite-depth code
%  which was PT's initial creation of the analytic continuation.
%  If possible this code should leverage some of the functions found there
%  instead of re-writing new code
%
%  Based on the Yyanis equations, we have 
%  
%           tau = I1 - I2 
%
%  YJL then explains later that I1 can be calculated as the sum of two
%  integrands, of which the second is the finite-depth integral. The first
%  integral he rewritees as I1 = T + Hilbert 
%
%  In later code the integrals are grouped like this: 
%
%  tau = H + (T - I2) 
%
%  In other words, we modify infinite depth codes only with the (T - I2)
%  integral, which is the "new" contribution. 

% So to make this code more accurate to the analytic continuation scheme,
% you should use the PVIGAUSS = 1 option. However, it is slightly less
% accurate than a midpoint scheme while on the axis it appears
USEPVIGAUSS = 0;

r0 = wave.r0;
N = length(phi);

% The solution vector consists of values of q and theta, and then mu and B
q = UU(1:N); tau = log(q);
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
q_m = D*q(:); q_m = q_m.'; tau_m = log(q_m);
theta_m = D*theta(:); theta_m = theta_m.';

% The singular equations we can evaluate using a midpoint scheme. However,
% note that we can override this with our PVI_GAUSS scheme, which seems to
% achieve exactly the same result; verify that bdintrhs_m is the the same
% as H
if USEPVIGAUSS ~= 1
    bdintrhs_m = 0*phi_m;
    for k=1:N-1
        f = theta.*cot(pi*(phi - phi_m(k)));
        g = 2*f;
        g(1) = g(1)/2;
        g(N) = g(N)/2;
        bdintrhs_m(k) = h/2*sum(g);    
    end
end

H_m = 0*phi_m;
I_m = 0*phi_m;

% Calculate the Hilbert + I-Integral

% For finite depth, we require an additional I integral given by (5.19) in
% YJL thesis
%
% Note that in the analytic continuation code, it is computed from the
% series solution, so we use this
%
% Iint = @(phi, f) qt_finite_bdint_Iint(phi, f, wave.an, wave.r0, wave.N);
%
% However, in this code, we just want to test the boundary integral, so we
% do it via the current mesh
Iint = @(f) tau - ...
  2*r0^2*(tau.*(r0^2 - cos(2*pi*(phi - f))) - theta.*sin(2*pi*(phi - f))) ...
              ./(r0^4 - 2*r0^2*cos(2*pi*(phi - f)) + 1);


for k = 1:length(phi_m)
    if USEPVIGAUSS == 1
        H_m(k) = qt_finite_getint_gauss(phi_m(k), wave);
    else
        % If you don't want to use the getint_gauss code, you can ge tthe
        % direct value here by commenting in the lines above 
        H_m(k) = bdintrhs_m(k);
    end
    if r0 == 0
        % Infinite depth has no integral
        I_m(k) = 0;
    else
        %I_m(k) = integral(@(phi)Iint(phi, phi_m(k)), -0.5, 0.5);        
        % I_m(k) = trapz(phi, Iint(phi, phi_m(k)));
        I_m(k) = trapz(phi, Iint(phi_m(k)));
    end
end
tau_m_from_integral = H_m + I_m;


%% Construct the nonlinear system 

F = zeros(2*N, 1);
% The first N-1 terms on the residual is from the integral
F(1:N-1) = log(q_m) - tau_m_from_integral;

%first the integral term.
m=(N+1)/2;
[~, Y] = qt_getxy(phi, q, theta);
Y_m = D*Y(:); Y_m = Y_m.';

% Bernoulli's
% Need derivative
% dqdphi_m = central_diff(q_m, phi_m);
% F(N:2*N-2) = q_m.^2.*dqdphi_m + 2*pi/mu*sin(theta_m);
%
% Note that I am using the integrated version of Bernoulli's. If I use
% the differentiated version, I don't end up using B. 
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

