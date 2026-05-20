function path = get_path(name)
    % Design standard paths
    % d is a list of parameters
    %
    % Use NaN to separate paths between upper/lower planes. The initial
    % condition for the contours will jump the NaN, so make sure pm (small
    % number) is chosen wisely (too close = numerical issues, too far = not
    % matched)
    %
    pm = 1i*1e-5;
    path.f = [];
    % This is an optional parameter that can be set for singularity location
    % If no value is set, set to Inf
    path.singstart = Inf;
    
    
    switch name
         case 'freesurf'
            phi0 = 0.4;
            f = [phi0 + pm, -phi0 + pm,-phi0 - pm, phi0 - pm]; 
            path.f = f;
        case 'boxaroundfree'
            phi0 = 0.4;
            psi0 = 1i*0.1;
            f = [phi0 + pm, phi0 + psi0, -phi0 + psi0, ...
                 -phi0 - psi0, phi0 - psi0, phi0 - pm];
            path.f = f;
        case 'box1pass'

            phi0 = 0.4;
            phi0right = 0.4;
            psi0 = 0.02*1i;
            f = [phi0 + pm, phi0 + psi0, psi0, -phi0 + psi0, -phi0 - psi0,...
                -psi0, phi0right - psi0, phi0right - pm, phi0right + pm, phi0right + psi0]; 
            path.f = f;
        case 'fluiddomain'
            phi0 = 0.4;
            psi0 = 0.05*1i;
            f = [phi0 + pm, phi0 - psi0, -phi0 - psi0, -phi0 -2*pm, phi0-2*pm];
            
        case 'box'
            phi0 = 0.4;
            psi0 = 0.2*1i;
            f = [phi0 + pm, phi0 + psi0, -phi0 + psi0, -phi0 + pm, NaN, ...
                -phi0 - pm, -phi0 - psi0, phi0 - psi0, phi0 - pm, NaN, ...
                phi0 + pm, phi0 + psi0, -phi0 + psi0, -phi0 + pm, NaN, ...
                -phi0 - pm, -phi0 - psi0, phi0 - psi0, phi0 - pm]; 
            path.f = f;
        case 'round.0.A'
            phi0 = 0.1;
            f = [phi0 + pm, phi0 + 1i, -phi0 + 1i, -phi0 + pm, ...
                phi0 + pm, phi0 + 1i, -phi0 + 1i, -phi0 + pm, ...
                phi0 + pm, phi0 + 1i, -phi0 + 1i, -phi0 + pm, ...
                phi0 + pm, phi0 + 1i, -phi0 + 1i, -phi0 + pm, ...
                phi0 + pm];
            path.f = f;
        case 'physicaldomain'
            phi0 = 0.3;
            psi0 = 0.1*1i;
            f = [phi0 + pm, phi0 - psi0, -phi0 - psi0, -phi0 + pm, phi0-pm]; 
            path.f = f;
        case 'testpath'
            phi0 = 0.4;
            phi1 = 0.3;
            phi2 = 0.2;
            phi3 = 0.1;
            psi0 = 0.1*1i;
            f = [phi0 + pm, phi0 + psi0, psi0, -phi0 + psi0, -phi0 - psi0,...
                -psi0, phi1 - psi0, phi1 + pm, phi1 + psi0, psi0, -phi1 + psi0, -phi1 - psi0,...
                -psi0, phi2 - psi0, phi2 + pm, phi2 + psi0, psi0, -phi2 + psi0, -phi2 - psi0,...
                -psi0, phi2 - psi0, phi2 + pm, phi3 + psi0, -phi0-psi0];
            path.f = f;
            
        case 'freesurf'
            phi0 = 0.3;
            psi0 = 0.0001*1i;
            f = [phi0 + pm, -phi0 + psi0, -phi0 + psi0, -phi0 - psi0, phi0-psi0, phi0-pm]; 
            path.f = f;
    end