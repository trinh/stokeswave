function [w, wp] = wfree_steep(f, bn, beta)
% WFREE constructs the series solution using coefficients from solseries

    
    t = exp(-2*pi*1i*f);
    A = 1 - beta*t;
    A13 = A.^(1/3);
    
    w = 1;
    for k = 1:length(bn)
        w = w + bn(k)*t.^k;
    end
    w = w.*A13;
    
    wp1 = 2*pi*1i*beta*t./(3*A).*w;    
    wp2 = 0;
    for k = 1:length(bn)
        wp2 = wp2 - 2*pi*1i*k*bn(k)*t.^k;
    end
    wp = wp1 + A13.*wp2;
    
            

end