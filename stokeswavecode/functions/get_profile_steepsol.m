function [phi, z, zp] = get_profile_steepsol(steepsolseries, n_phi)
    % GET_PROFILE_STEEPSOL calculates the z(phi) and dz/dphi values for a
    % steepsolseries structure. It uses a special gridding near the crest
    
    if nargin < 2
        n_phi = 81;
    end
    
    L = 0.5;        % Half-domain length
    pow = 1/3;      % Assume x' and y' ~ f^(1/3)
    
    mid = (n_phi + 1)/2;
    
    % Create a phi vector that is scaled so that f^(1/3) is equally spaced
    phi = nan*ones(1, n_phi);
    s = linspace(0, L.^pow, mid);
    phi(mid+1:end) = s(2:end).^(1/pow);
    phi(1:mid-1) = -fliplr(phi(mid+1:end));
    phi(mid) = 0;
    
    
    [w, wp] = vb_wfree_steep(phi, steepsolseries.bn, steepsolseries.beta);
    
    zp = 1./w;
    
    x_phi = real(zp);
    y_phi = imag(zp);
    
    % Integrate from zero to get (x, y)
    x = 0*w;
    y = 0*w;
    
    for j = mid+1:length(w)
        x(j) = x(j-1) + (phi(j) - phi(j-1))*(x_phi(j) + x_phi(j-1))/2;
        y(j) = y(j-1) + (phi(j) - phi(j-1))*(y_phi(j) + y_phi(j-1))/2;
    end
    for j = mid-1:-1:1
        x(j) = x(j+1) + (phi(j) - phi(j+1))*(x_phi(j) + x_phi(j+1))/2;
        y(j) = y(j+1) + (phi(j) - phi(j+1))*(y_phi(j) + y_phi(j+1))/2;
    end
    
    z = x + 1i*y;
end