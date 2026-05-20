function theta = get_theta_series(f, an, r0, N)
    % QT_SERIES calculates (q, theta) from series data
    %
    % 2026-05-20: I want to phase out this function. It should be embedded
    % into the VBStokesWave class so that we don't make mistakes
    
    error('This function should not be used; see notes')

    zp = 1;
    for k = 1:(N - 1)
        zp = zp + an(k)*(exp(-2i*pi*f*k) + r0^(2*k)*exp(-2i*pi*f*(-k)));
    end
    
    % From the series data need to get the speed and angle
    % q = 1./abs(zp);
    theta = atan(imag(zp)./real(zp));    
end