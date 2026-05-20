function [Cnew, s] = getComplexCurveByArcLength(C, N_interp, reimname, arcmeasure)
% PARAMARC Reparameterizes a curve by arclength
%
%   PHT NOTE: used in stock functions
%
%   Given a N-D curve via vectors C (2xm) we create a parameterization so
%   that C = C(alpha) where alpha begins at 0 and goes to
%   1, but is scaled according to the segment lengths in between
%
%   arcmeasure = 'flatten' means measure arclength only in the C(1,:) dimension
%
%   Then we interpolate using N_interp points.
%
%   The arclength uses the 3D point (x, y, reimname(f))

reim = str2func(reimname);

[N_dim, N_length] = size(C);
s = linspace(0, 1, N_interp);

if N_dim ~= 2
    keyboard
end

if nargin == 3
    arcmeasure = '3d';
end

% Form 3d contour to measure the arclength
switch arcmeasure
    case 'flatten'
        C_arc = [real(C(1,:)); imag(C(1,:)); 0*reim(C(2,:))];
    otherwise
        C_arc = [real(C(1,:)); imag(C(1,:)); reim(C(2,:))];
end

%update arclength
alpha = zeros(1, N_length);
alpha(1)=0;
for j = 2:N_length
    alpha(j) = norm(C_arc(:,j)-C_arc(:,j-1));
end

alpha = cumsum(alpha);
alpha = alpha/alpha(end); %normalized arclength parametrization

% Anytime there is a zero difference, the next point repeats, so remove
dalpha = diff(alpha);
ind = find(dalpha == 0);
alpha(ind+1) = [];
C(:,ind+1) = [];

%reparameterization
Cnew = nan(N_dim, N_interp);
for j = 1:N_dim
    Cnew(j,:) = interp1(alpha, C(j,:), s,'spline');
end

% Test
%{
figure
plot3(C_arc(1,:), C_arc(2,:), C_arc(3,:), 'bo-');
hold on
plot3(real(Cnew(1,:)), imag(Cnew(1,:)), reim(Cnew(2,:)), 'r.-');
hold off
keyboard
%}

end

