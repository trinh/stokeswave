function theta = get_theta_series(f, solseries)
% QT_SERIES calculates (q, theta) from series data

    an = solseries.an;
    zp = 1;
    for k = 1:(solseries.N - 1)
        zp = zp + an(k)*exp(-2*pi*k*1i*f);
    end 
            
    % From the series data need to get the speed and angle
    % q = 1./abs(zp);
    theta = atan(imag(zp)./real(zp));

end