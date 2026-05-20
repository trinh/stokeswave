function q = get_q_series(f, solseries)
% QT_SERIES calculates (q, theta) from series data

    an = solseries.an;
    dzdf = 1;
    for k = 1:(solseries.N - 1)
        dzdf = dzdf + an(k)*exp(-2*pi*k*1i*f);
    end 
            
    % From the series data need to get the speed and angle
    q = 1./abs(dzdf);
    % theta = atan(imag(zp)./real(zp));

end