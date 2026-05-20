function [z, zp] = zfree(f, solseries)
% ZFREE constructs the series solution using coefficients from solseries

    an = solseries.an;
    z = f;
    zp = 1;
    for k = 1:(solseries.N - 1)
        z = z + ((1i*an(k))/(2*pi*k))*exp(-2*pi*k*1i*f);
        zp = zp + an(k)*exp(-2*pi*k*1i*f);
    end
    z = z - 1i*solseries.ycrest;
            

end