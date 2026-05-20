function zp = zpfree(f, solseries)
% ZPFREE constructs the series solution using coefficients from solseries
% 
% Written so we can easily get zp as one-output

    an = solseries.an;
    zp = 1;
    for k = 1:(solseries.N - 1)
        zp = zp + an(k)*exp(-2*pi*k*1i*f);
    end 
            

end