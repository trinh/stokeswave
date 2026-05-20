function F = series_bern_steep2(Y, ep, N, usebeta)

% a the vector of unknowns. I.e. a_n coefficients and mu and beta.
% ep given wave "steepness"
% N step size

% b = Y(1:N-1);
% mu = Y(N);
% % beta = abs(Y(N+1));
% beta = 0;

if usebeta == 0
    b = Y(1:N);
    mu = Y(N+1);
    beta = 0;
elseif usebeta == 1
    b = Y(1:N-1);
    beta = Y(N);
    mu = Y(N+1);
elseif usebeta == 2
    b = Y(1:N);
    mu = Y(N+1);
    beta = 1;
end
% beta = abs(beta);

F = zeros(1, N+1);
phi = (2*(1:N)-1)/(4*N);
[w, wp] = wfree_steep(phi, b, beta);
% wp = central_diff(w, phi);
u = real(w); up = real(wp);
v = -imag(w); vp = -imag(wp);
Wabs = sqrt(u.^2 + v.^2);
WabsWabsp = up.*u + vp.*v;
F(1:N) = WabsWabsp - 2*pi/mu*imag(w)./Wabs.^2;

[w0, ~] = wfree_steep(0, b, beta);
[wend, ~] = wfree_steep(0.5, b, beta);
F(N+1) = 1 - (abs(w0)^2)*(abs(wend)^2) - ep^2;

%  plot(abs(F))
% % keyboard
%  drawnow
