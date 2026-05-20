function path = get_path_ep0p3(name)
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
    
    fA = 1i*0.176;
    fB = 0.79 - 0.55i;
    fBl = -real(fB) + 1i*imag(fB);
    fC = 1.257 - 0.837i;
    fCl = -real(fC) + 1i*imag(fC);
    % fD = 1.14 - 0.8i
    
    xmax = 1.8;
    phi0 = 0.1;
    psi0 = 1i*(imag(fA) + min(2*imag(fA), 0.1));
    psifudge = 0.01i;
    
    circle_points = 30;
    circle_tours = 1;
    circle_rad = 0.1;
    
    % One tour
    tour = [phi0 + pm, phi0 + psi0, -phi0 + psi0, -phi0 + pm, NaN, ...
        -phi0 - pm, -phi0 - psi0, phi0 - psi0, phi0 - pm];
    
    % Tour around B
    if real(fB) > 1
        warning('the B-tour will also include A copies')
    end
    phi0B = real(fB) + circle_rad;
    psi0B = 1i*(imag(fB) + circle_rad);
    tourB = [0 + pm, phi0B + pm, phi0B + psi0B, -phi0B + psi0B, -phi0B + pm, ...
            NaN, -phi0B - pm, -phi0B - psi0B, phi0B - psi0B, phi0B - pm];
    
    switch name
        
        % ========================================================
        %                  DUAL MESHES
        %
        %               Used for publications
        % * It seems we need at least two non-trivial elements
        %    following a NaN declaration!
        % ========================================================
        case 'dualmesh.0'
            base = [phi0 + pm, -xmax + pm];
            path.up = [base, 999];
            path.down = [base, NaN, -xmax - pm, -xmax - psifudge, -xmax - pm, 999];
            path.f = path.up(1:end-1);
        case 'dualmesh.A'
            base = [tour(1:4), -xmax + pm];
            path.up = [base, 999];
            path.down = [base, NaN, -xmax - pm, -xmax - psifudge, -xmax - pm, 999];
            path.f = path.up(1:end-1);
        case 'dualmesh.A1'
            base = [0 + pm, 1+tour(1:4), -xmax + pm];
            path.up = [base, 999];
            path.down = [base, NaN, -xmax - pm, -xmax - psifudge, -xmax - pm, 999];
            path.f = path.up(1:end-1);
        case 'dualmesh.A1A1_'
            base = [0 + pm, 1+tour, -xmax - pm];
            path.down = [base, 999];
            path.up = [base, NaN, -xmax + pm, -xmax + psifudge, -xmax + pm, 999];
            path.f = path.down(1:end-1);
        case 'dualmesh.A1A_'
            base = [0 + pm, 1+tour(1:4), NaN, 1 - phi0 - pm, tour(6:end), -xmax - pm];
            path.down = [base, 999];
            path.up = [base, NaN, -xmax + pm, -xmax + psifudge, -xmax + pm, 999];
            path.f = path.down(1:end-1);            
        case 'dualmesh.An1'
            base = [0 + pm, -1+tour(1:4), -xmax + pm];
            path.up = [base, 999];
            path.down = [base, NaN, -xmax - pm, -xmax - psifudge, -xmax - pm, 999];
            path.f = path.up(1:end-1);
        case 'dualmesh.AA_'
            base = [tour, -xmax - pm];
            path.down = [base, 999];
            path.up = [base, NaN, -xmax + pm, -xmax + psifudge, -xmax + pm, 999];
            path.f = path.down(1:end-1);
        case 'dualmesh.AA_B_'
            base = [tour(1:8), real(fB) - circle_rad + 1i*imag(fB)];
            circtour = gencirc(fB, base(end), circle_points, circle_tours);
            base = [base, circtour, real(base(end)) - pm, -xmax - pm];
            path.down = [base, 999];
            path.up = [base, NaN, -xmax + pm, -xmax + 0.05i, -xmax + pm, 999];
            path.f = path.down(1:end-1);
        case 'dualmesh.AA_B_^'
            base = [tour, real(fBl) - pm]; 
            base = [base, fBl + 1i*circle_rad];
            circtour = gencirc(fBl, base(end), circle_points, circle_tours);
            base = [base, circtour, real(base(end)) - pm, -xmax - pm];
            path.down = [base, 999];
            path.up = [base, NaN, -xmax + pm, -xmax + 0.05i, -xmax + pm, 999];
            path.f = path.down(1:end-1);
        case 'dualmesh.AA_B_^B_'
            base = [tourB, -xmax - pm];
            path.down = [base, 999];
            path.up = [base, NaN, -xmax + pm, -xmax + 0.05i, -xmax + pm, 999];
            path.f = path.down(1:end-1);
        case 'dualmesh.AA_A1'
            base = [tour, NaN, phi0 + pm, 1 + tour(1:4), -xmax + pm]; 
            path.down = [base, 999];
            path.up = [base, NaN, -xmax + pm, -xmax + psifudge, -xmax + pm, 999];
            path.f = path.down(1:end-1);
        case 'dualmesh.AA_AA_'
            base = [tour, NaN, tour, -xmax - pm];
            path.down = [base, 999];
            path.up = [base, NaN, -xmax + pm, -xmax + psifudge, -xmax + pm, 999];
            path.f = path.down(1:end-1);
        case 'dualmesh.AA_AA_B_'
           base = [tour, NaN, tour(1:8), real(fB) - circle_rad + 1i*imag(fB)];
            circtour = gencirc(fB, base(end), circle_points, circle_tours);
            base = [base, circtour, real(base(end)) - pm, -xmax - pm];
            path.down = [base, 999];
            path.up = [base, NaN, -xmax + pm, -xmax + 0.05i, -xmax + pm, 999];
            path.f = path.down(1:end-1);
        case 'dualmesh.AA_AA_B_^'
            base = [tour, NaN, tour, real(fBl) - pm]; 
            base = [base, fBl + 1i*circle_rad];
            circtour = gencirc(fBl, base(end), circle_points, circle_tours);
            base = [base, circtour, real(base(end)) - pm, -xmax - pm];
            path.down = [base, 999];
            path.up = [base, NaN, -xmax + pm, -xmax + 0.05i, -xmax + pm, 999];
            path.f = path.down(1:end-1);
        case 'dualmesh.AA_AA_C_'
            base = [tour, NaN, tour, NaN, tour(1:7)];
            base = [base, real(fC) + 1i*imag(base(end)), fC + 1i*circle_rad];
            circtour = gencirc(fC, base(end), circle_points, circle_tours);
            base = [base, circtour, real(base(end)) - pm, -xmax - pm];
            path.down = [base, 999];
            path.up = [base, NaN, -xmax + pm, -xmax + 0.05i, -xmax + pm, 999];
            path.f = path.down(1:end-1);            
        case 'dualmesh.AA_AA_C_^'
            base = [tour, NaN, tour, NaN, tour, real(fCl) - pm];
            base = [base, fCl + 1i*circle_rad];
            circtour = gencirc(fCl, base(end), circle_points, circle_tours);
            base = [base, circtour, real(base(end)) - pm, -xmax - pm];
            path.down = [base, 999];
            path.up = [base, NaN, -xmax + pm, -xmax + 0.05i, -xmax + pm, 999];
            path.f = path.down(1:end-1);
        case 'dualmesh.AA_AA_AA_'
            base = [tour, NaN, tour, NaN, tour, -xmax - pm];
            path.down = [base, 999];
            path.up = [base, NaN, -xmax + pm, -xmax + psifudge, -xmax + pm, 999];
            path.f = path.down(1:end-1);
            
    end
end % end main function

function z = gencirc(zcentre, z1, numpts, numtours)
    % GENCIRC generates points around a circle
    rad = abs(z1 - zcentre);
    theta0 = angle(z1 - zcentre);
    t = linspace(0, numtours*2*pi, numtours*numpts + 1);
    % Remove the first point, which repeats
    t = t(2:end);
    z = zcentre + rad.*exp(1i*(theta0 + t));
end
function z = genspiral(z0, z1, numpts, numtours, fac)
    % GENSPIRAL generates points around a spiral
    rad = abs(z1 - z0);
    theta0 = angle(z1 - z0);
    t = linspace(0, numtours*2*pi, numtours*numpts + 1);
    t = t(2:end);
    z = z0 + (rad + fac*rad*sin(t/numtours)).*exp(1i*(theta0 + t));
end
function z = genellipse(z1, zs, a, anglerange, numpts)
    % GENELLIPSE creates an ellipse where zs is the centre and z1 is the
    % point on the b-semi-major axis. The other axis (a) is created
    % perpendicular to the vector joining zs and z1
    bv = z1 - zs;
    u = 1i*bv; u = u/abs(u);
    
    t = linspace(anglerange(1), anglerange(2), numpts + 1);
    t = t(2:end);
    
    z = zs + a*u*sin(t) + bv*cos(t);
end
function z = connectcircle(z0, zs, numpts, numtours)
    % CONNECTCIRCLE connects the point (x0, 0) to a circle
    
    R = imag(z0);   % maximal distance
    v = (z0 - zs);  % vector from centre to starting point
    v = v/abs(v);   % normalize
    
    z1 = zc + (R/2)*v;
    
    znew = gencirc(zs, z1, numpts, numtours);
    z = [z1, znew];
end