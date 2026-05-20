function [gp, p] = gbern_lower2upper(f, g, solseries)
    N = solseries.N;
    B = solseries.B;
    mu = solseries.mu;
    an = solseries.an;
    
    p = ones(size(f));
    gp = 0*f;
    for k = 1:N-1
        p = p + an(k).*exp(2*pi*1i*k*f);
    end
    
    gp = mu./(2*p.*(mu*B + 1i*pi*g)) - p;    
    
    
end