function F = bdint_path_bdfunc1(w, Y, solbdint, solseries, meshdatbd)
% BDFUNC_COMPLEX is the complex version to find the analytic continuation
% 

    % always set a = 1 for the initial integration upwards
    a = 1;
    
    sheet_number = 1;
    F = zeros(2*sheet_number, 1);
    x = Y(1:sheet_number);
    y = Y(sheet_number+1:2*sheet_number);
    
    
    % Technically, we should reverse a = -1 for initial step into the LHP
    % However, the inLHP flag should handle this
    inLHP = 0*w;
    inLHP(imag(w) < 0) = 1;        
    
    % Must assume that input is of size (2, N_w)
    x = Y(1, :);
    y = Y(2, :);
    
    B = solbdint.B;
    mu = solbdint.mu;
    
    % G = getint(w, solbdint);
    G = getint_gauss(w, solseries);
    if abs(imag(w)) < 1e-10
        keyboard
    end
    
    % How should we get beta in the LHP?
    betafree = 0*w;
    for j = 1:length(w)
    %     [~, zp] = zfree(w, solseries);
    %     betafree = imag(zp);
        if inLHP(j) == 1
            betafree(j) = interp2(real(meshdatbd.wmat), imag(meshdatbd.wmat), meshdatbd.ypmat, real(w(j)), imag(w(j)));
%             betafree(j) = interp2(real(meshdatbd.wmat), imag(meshdatbd.wmat), imag(meshdatbd.xpmat + 1i*meshdatbd.ypmat), real(w(j)), imag(w(j)));
%             [~, zp] = zfree(w, solseries);
%             betafree = imag(zp);
%             betafree(j) - betafree2
%             keyboard
        else
            betafree(j) = 0;
        end
    end
    
    
    yp = a*1i./(2*(G - 2i*betafree*inLHP)).* ...
            ((G - 2i*betafree*inLHP).^2 - (1/2)./(B - 2*pi/mu*y));
    
    xp = G + a*1i*yp + inLHP*(-2i*betafree);
%     keyboard
    F(1,:) = xp;
    F(2,:) = yp;    
    
end