function [X, Y] = qt_getxy(phi, q, theta)
% QT_GETXY gets the X and Y values from normal integration from (q, theta)
%
% Assumes that phi is odd and zero is in the middle

N = length(phi);
m = (N+1)/2;

Y = 0*phi;
h = phi(2) - phi(1);
beta = sin(theta)./q;
for j = m:N-1
    Y(j+1) = Y(j) + h/2*(beta(j+1) + beta(j));
end
for j = m:-1:2
    Y(j-1) = Y(j) - h/2*(beta(j-1) + beta(j));
end

X = cumtrapz(phi, cos(theta)./q) - 0.5;