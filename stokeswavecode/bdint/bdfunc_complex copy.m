function F = bdfunc_complex(w, Y, solbdint, solseries)
% BDFUNC_COMPLEX is the complex version to find the analytic continuation

    % dir = +1 for upwards, -1 for downwards
    % always set a = 1 for the initial integration upwards
    a = 1;
    if imag(w) < 0
        inLHP = 1;
    else
        inLHP = 0;
    end

    x = Y(1);
    y = Y(2);
    
    w = w(:).';
    
    B = solbdint.B;
    mu = solbdint.mu;
    
    G = getint(w, solbdint);
    
    [~, zp] = zfree(w, solseries);
    betafree = imag(zp);
    yp = a*1i./(2*(G - 2i*betafree*inLHP)).* ...
            ((G - 2i*betafree*inLHP).^2 - (1/2)./(B - 2*pi/mu*y));
    
    xp = G + a*1i*yp + inLHP*(-2i*betafree);
    F = [xp; yp];
    
end