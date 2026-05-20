function zpp = zppfree(f, solseries)
% ZPPFREE constructs the series solution using coefficients from solseries
% 
% Written so we can easily get zpp as one-output

    an = solseries.an;
    zpp = 0;
    for k = 1:(solseries.N - 1)
        zpp = zpp + an(k)*(-2i*pi*k)*exp(-2i*pi*k*f);
    end 
            

end