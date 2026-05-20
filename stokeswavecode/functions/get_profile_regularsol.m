function [phi, z, zp] = get_profile_regularsol(solseries, n_phi)
    % GET_PROFILE_REGULARSOL calculates the z(phi) and dz/dphi values for a
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
    
    [z, zp] = zfree(phi, solseries);    
    
end