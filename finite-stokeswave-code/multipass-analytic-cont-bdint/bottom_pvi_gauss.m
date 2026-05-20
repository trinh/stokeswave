function [I, collapse] = pvi_gauss(f, x0, a, b, r0)
% PVI_GAUSS calculates a singular integral in two ways:
%   When the singularity is away from the axis, it simply integrates using
%   the normal adaptive integration scheme
% 
%   When the singularity is near or on the axis, it first pushes the
%   singularity right down to the axis, and then does the scheme by
%   splitting the integral into three sections. The first/third two sections are
%   well behaved. The middle is via 8-point Gauss-Legendre
%
%   2016-03-25: Edit to allow specification of tolerances
%   2015-08-01: (Approx.) Coding for the Stokes wave

abstol      = 1e-10;    % Integral abs tol
reltol      = 1e-6;     % Integral rel tol
collapsetol = 1e-8;     % Collapse tolerance
delta       = 1e-4;     % Distance for the residue


myint = @(a1, a2)integral(f, a1, a2, 'AbsTol', abstol, 'RelTol', reltol);
Q = -log(r0)/(2*pi);
% Anything that lies within mytol of the axis should be collapsed, and then
% we manually add in the residue

if abs(imag(x0)+2*Q) < collapsetol
    x0 = real(x0)-2*Q*1i;
    collapse = 1;
else
    collapse = 0;
end

if (imag(x0) == -2*Q*1i) && (real(x0) > a) && (real(x0) < b)
    
    % Integrate in three parts
    % Take default value of delta 
    % I'm not sure how to choose this effectively    

    % First part from a to x0 - L
    I1 = myint(a, real(x0) - delta);

    % Second part from x0 - L to x0 + L
    % Gauss-Legendre 8 - point
    s = [0.960289856497536 ...
       0.796666477413627 ...
       0.525532409916329 ...
       0.183434642495650 ...
      -0.183434642495650 ...
      -0.525532409916329 ...
      -0.796666477413627 ...
      -0.960289856497536];
    wt = [0.101228536290377 ...
       0.222381034453374 ...
       0.313706645877887 ...
       0.362683783378362 ...
       0.362683783378362 ...
       0.313706645877887 ...
       0.222381034453374 ...
       0.101228536290377];
    I2 = delta*sum(wt.*f(delta*s + x0));
    
    % [x, w] = lgwt(8, x0 - delta, x0 + delta);
    % I2 = sum(w.*f(x));    

    % Third part from x0 + L to b
    I3 = myint(real(x0) + delta, b);

    I = I1 + I2 + I3;
else
    I = myint(a, b);    
end