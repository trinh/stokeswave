 function Iintegrand = qt_bdint_Iint(phi, f, an, r0, N)
    % Computes the secondary integrand (tau - r0^2*[...])
    %   Note that phi is the integration variable and f is a parameter
    %
    %   2016-12-01: Coding for YJL
    
    zp = 1;
    for k = 1:N - 1
        zp = zp + an(k)*(exp(-2i*pi*phi*k) + r0^(2*k)*exp(-2i*pi*phi*(-k)));
    end
    
    % From the series data need to get the speed and angle
    tau = log(1./abs(zp));
    theta = atan(imag(zp)./real(zp));    
    
    Iintegrand = tau;
end