function path = get_path_2026_05_08(name)
    % Design standard paths
    % d is a list of parameters
    %
    % Use NaN to separate paths between upper/lower planes. The initial
    % condition for the contours will jump the NaN, so make sure pm (small
    % number) is chosen wisely (too close = numerical issues, too far = not
    % matched)
    %
    %     Notation:
    %         round.[A].[B] gives a path round singularity [B] whose final
    %         Riemann sheet is [A]
    %
    %         dualmesh.[A] gives a dual mesh structure (only curently used in
    %         qt_meshsurface_phi.m
    %
    %         bang.[A].[B] gives path that hopes to directly hit (and so ODE
    %         fails) singularity [B] and ends up on final Riemann sheet [A]
    %
    %         circle.[A].[B] is like round but uses a precise circle around a
    %         known singularity
    %
    %    Notes:
    %         Use periods (.) to add notes to end of the above notation.
    %
    
    pm = 1i*1e-5;
    path.f = [];
    % This is an optional parameter that can be set for singularity location
    % If no value is set, set to Inf
    path.singstart = Inf;
    
    
    switch name
        case 'bang.0.A'
            f = [0, 1i];
            path.f = f;
            % ========================================================
            %                  ROUND PATHS
            %
            %    These are straight-line paths
            %
            % ========================================================
            
            % a = A, A = \uline{A}, B^ = Bl
        case 'round.0.A'
            phi0 = 0.1;
            f = [phi0 + pm, phi0 + 1i, -phi0 + 1i, -phi0 + pm, ...
                phi0 + pm, phi0 + 1i, -phi0 + 1i, -phi0 + pm, ...
                phi0 + pm, phi0 + 1i, -phi0 + 1i, -phi0 + pm, ...
                phi0 + pm, phi0 + 1i, -phi0 + 1i, -phi0 + pm, ...
                phi0 + pm];
            path.f = f;
    end
    