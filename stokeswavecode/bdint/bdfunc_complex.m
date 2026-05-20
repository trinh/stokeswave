function F = bdfunc_complex(w, Y, solbdint, solseries, a, inLHP)
% BDFUNC_COMPLEX is the complex version to find the analytic continuation
% 
%   To handle vector input: w is size 1 x N 
%   F should be size (2, N)

    if size(Y, 1) ~= 2 || size(w, 1) ~= 1
        warning('Y should be 2 x length(w)')
        keyboard
    end

    N_w = length(w);
    F = zeros(2, N_w);
    
    % always set a = 1 for the initial integration upwards
    if nargin < 5
        a = 1;
    end
    
    % Technically, we should reverse a = -1 for initial step into the LHP
    % However, the inLHP flag should handle this
    if nargin < 6
        inLHP = 0*w;
        inLHP(imag(w) < 0) = 1;        
    end

    % Must assume that input is of size (2, N_w)
    x = Y(1, :);
    y = Y(2, :);
    
    B = solbdint.B;
    mu = solbdint.mu;
    
    % G = getint(w, solbdint);  
    G = 0*w;
    for j = 1:length(w)
        G(j) = getint_gauss(w(j), solseries);
    end
    
    [~, zp] = zfree(w, solseries);
    betafree = imag(zp);
    yp = a*1i./(2*(G - 2i*betafree*inLHP)).* ...
            ((G - 2i*betafree*inLHP).^2 - (1/2)./(B - 2*pi/mu*y));
    
    xp = G + a*1i*yp + inLHP*(-2i*betafree);
    F(1,:) = xp;
    F(2,:) = yp;    
    
end