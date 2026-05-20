function path = get_path_finite(name)
    %
    % WARNING: THIS IS AN INFERIOR VERSION OF A TEMPLATE get_path THAT PHIL
    % WROTE BEFORE. THIS WILL LIKELY BE ABANDONED 2026 MAY TO ALLOW
    % GET_PATH TO BE SPECIFIED WITHIN INDIVIDUAL STUDY FOLDERS
    %
    %
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
        case 'box1pass'
            phi0 = 0.1;
            psi0 = 0.2*1i;
            f = [phi0 + pm, phi0 + psi0, -phi0 + psi0, -phi0 + pm, NaN, ...
                -phi0 - pm, -phi0 - psi0, phi0 - psi0, phi0 - pm]; 
            path.f = f;
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
            f = [phi0 - pm, phi0 - psi0, -phi0 - psi0, -phi0 + pm, phi0-pm]; 
            path.f = f;
    end