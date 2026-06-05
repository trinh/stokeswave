function out = getIBad(xi,Y,h,Z,L)
% A function to compute the integral I
% Input: xi = 
%           L = truncation limits of infinite integral

Mk=1i*1024*pi/0.5;
Yhat = fft(Y); 
Yxi = real(ifft(Mk.*Yhat));

    [ny, nx] = size(Z); 
    out = zeros(ny,nx);
    for county = 1:ny
        for countx = 1:nx
            ZZ = Z(county,countx);
            YxiCon = @(sig) getYxiCon(xi,Yxi,sig);
            I = @(sig) YxiCon(sig)/(2*h).*coth(pi*(sig-ZZ)/(2*h));
            % Itaylor = @(sig) YxiCon(sig)/pi;
            out(county, countx) = pvi_gauss(I,ZZ,-L,L);
            % out(county,countx) = pvi_gauss2(I,ZZ,-L,L, Itaylor);
        end
   end
end

function out = getYxiCon(xi, Yxi, sig)
    % getYxiCon
    % Interpolates the value of Y(xi) at a given angular parameter sig,
    % accounting for periodicity on the interval [-0.5, 0.5).

    % Shift sig into the interval [-0.5, 0.5) using modular arithmetic.
    % This ensures periodic behavior without discontinuities.
    sig_shift = mod(sig + 0.5, 1) - 0.5;

    % Because xi and Yxi represent periodic data, we manually append the first
    % value at xi = -0.5 to the end at xi = 0.5. This avoids interpolation gaps.
    % Append Yxi(-0.5) to match the periodic boundary.
    Yxi = [Yxi, Yxi(1)];
    
    % Append the corresponding xi value defining the full periodic span.
    xi = [xi, 0.5];

    % Interpolate Y(xi) at the shifted sig location.
    out = interp1(xi, Yxi, sig_shift);

    % If interpolation fails (should not happen if data is periodic and monotonic),
    % drop to debugging mode.
    if isnan(out)
        keyboard
    end
end