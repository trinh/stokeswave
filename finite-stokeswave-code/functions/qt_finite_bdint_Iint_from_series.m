function Iintegrand = qt_finite_bdint_Iint_from_series(phi, f, an, r0, N)
    % Computes the secondary integrand (tau - r0^2*[...])
    %   Note that phi is the integration variable and f is a parameter
    %
    %   2016-12-01: Coding for YJL
    %   2026-05-19: Renaming finite so as to avoid catastrophic problems
    %
    % --------------------------------------------------------------------
    %  Note that this file was renamed to emphasise that it is
    %  series-dependent; this allows adaptive integration 

    zp = 1;
    for k = 1:N - 1
        zp = zp + an(k)*(exp(-2i*pi*phi*k) + r0^(2*k)*exp(-2i*pi*phi*(-k)));
    end

    % From the series data need to get the speed and angle
    tau = log(1./abs(zp));
    theta = atan(imag(zp)./real(zp));

    % Correction 2026 May 19: I think we need factor of 2
    Iintegrand = tau - ...
        2*r0^2*(tau.*(r0^2 - cos(2*pi*(phi - f))) ...
        - theta.*sin(2*pi*(phi - f))) ...
        ./(r0^4 - 2*r0^2*cos(2*pi*(phi - f)) + 1);
end