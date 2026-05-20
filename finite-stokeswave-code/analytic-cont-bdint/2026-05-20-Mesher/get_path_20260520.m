function path = get_path_20260520(name)
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
        % ========================================================
        %                   BANG PATHS
        %
        % Paths that should bang into the singularities so as to provide a
        % crude check of the singularity values
        %   Notation - bang.[sheet].[singularity]
        % ========================================================
        case 'bang.0.A'
            f = [0, 1i];
        case 'bang.A.A_'
            f = [0.2, 0.2 + 1i, -0.2 + 1i, -0.2 + pm, 0 + pm, NaN, 0 - pm, -0.3i];
            
            % ========================================================
            %                   MESHING PATHS
            %
            % Specify paths to mesh the plane.
            % Use 999 to indicate where plane is first meshed
            % The coordinate previous to 999 will be the start point
            % Mesh in UHP/LHP depending on last coordinate
            %
            % ========================================================
            
            % Used to mesh the first UHP
            %         case 'mesh_n1_up'
            %             edgesize = 1.5;
            %             f = [-0.4 + pm, -0.4 + 1i, 0.4 + 1i, 0.4 + pm, -edgesize + pm, 999];
            %         case 'mesh_n1_up_test'
            %             f = [-0.4 + pm, -0.4 + 1i, 0.4 + 1i, 0.4 + pm, -0.4 + pm, ...
            %                 -0.4 + 1i, 0.4 + 1i];
            %             % Used to mesh the first unphysical fluid
            %         case 'mesh_n1_down'
            %             edgesize = 1.5;
            %             f = [-0.4 + pm, -0.4 + 1i, 0.4 + 1i, 0.4 + pm, NaN, ...
            %                 0.4 - pm, -edgesize - pm, 999];
            %         case 'mesh_n1_down_test'
            %             f = [-0.4 + pm, -0.4 + 1i, 0.4 + 1i, 0.4 + pm, NaN, ...
            %                 0.4 - pm, 0.4 - 1i, -0.4 - 1i];
            %         case 'mesh_n1n1_down'
            %             edgesize = 1.5;
            %             f = [-0.4 + pm, -0.4 + 1i, 0.4 + 1i, 0.4 + pm, NaN, ...
            %                 0.4 - pm, 0.4 - 1i, ...
            %                 -0.4 - 1i, -0.4 - pm, -edgesize - pm, 999];
            %         case 'mesh_n1n1_down_test'
            %             f = [-0.4 + pm, -0.4 + 1i, 0.4 + 1i, 0.4 + pm, NaN, ...
            %                 0.4 - pm, 0.4 - 1i, ...
            %                 -0.4 - 1i, -0.4 - pm, -1.5 - pm, -1.5 - 1i, 1.5 - 1i];
            %         case 'mesh_n1n1_up'
            %             edgesize = 1.5;
            %             f = [-0.4 + pm, -0.4 + 1i, 0.4 + 1i, 0.4 + pm, NaN, ...
            %                 0.4 - pm, 0.4 - 1i, ...
            %                 -0.4 - 1i, -0.4 - pm, NaN, -0.4 + pm, ...
            %                 -edgesize + pm, 999];
            
        case 'mesh.0.up'
            xmax=0.8;
            f = [-xmax + pm, xmax + pm, 999];
            path.f = f;
        case 'mesh.A.up'
            xmax=0.8;
            f = [-xmax + pm, -xmax + 1i, xmax + 1i, xmax + pm, -xmax + pm, 999];
        case 'mesh.A.down'
            xmax=0.8;
            f = [-xmax + pm, -xmax + 1i, xmax + 1i, xmax + pm, NaN, ...
                xmax - pm, -xmax - pm, 999];
            path.f = f;
        case 'mesh_BE_down'
            xmax=0.9;
            f = [-xmax + pm, -xmax + 1i, xmax + 1i, xmax + pm, NaN, ...
                xmax - pm, xmax - 1i, -0.1 - 1i, -0.1 - pm, -xmax - pm, 999 ];
        case 'mesh_BE_up'
            xmax=0.9;
            f = [-xmax + pm, -xmax + 1i, xmax + 1i, xmax + pm, NaN, ...
                xmax - pm, xmax - 1i, -0.1 - 1i, -0.1 - pm, NaN, ...
                -0.1 + pm, -xmax + pm, 999];
        case 'mesh_BED_down'
            xmax=0.9;
            f = [-xmax + pm, -xmax + 1i, xmax + 1i, xmax + pm, NaN, ...
                xmax - pm, xmax - 1i, -xmax - 1i, -xmax - pm, 999];
        case 'mesh_BED_up'
            xmax=0.9;
            f = [-xmax + pm, -xmax + 1i, xmax + 1i, xmax + pm, NaN, ...
                xmax - pm, xmax - 1i, -xmax - 1i, -xmax - pm, NaN, ...
                -xmax + pm, -xmax + 1i, -xmax + pm, 999];
        case 'mesh_BEDB_up'
            xmax=0.8;
            f = [-xmax + pm, -xmax + 1i, xmax + 1i, xmax + pm, NaN, ...
                xmax - pm, xmax - 1i, -xmax - 1i, -xmax - pm, NaN, ...
                -xmax + pm, -xmax + 1i, 0.1 + 1i, 0.1 + pm, -xmax + pm, 999];      %
        case 'mesh_BEDB_down'
            xmax=0.8;
            f = [-xmax + pm, -xmax + 1i, xmax + 1i, xmax + pm, NaN, ...
                xmax - pm, xmax - 1i, -xmax - 1i, -xmax - pm, NaN, ...
                -xmax + pm, -xmax + 1i, 0.1 + 1i, 0.1 + pm, NaN, ...
                0.1 - pm, -xmax - pm, 999];
            
        case 'mesh_BEDE_down'
            xmax=0.8;
            f = [-xmax + pm, -xmax + 1i, xmax + 1i, xmax + pm, NaN, ...
                xmax - pm, xmax - 1i, -xmax - 1i, -xmax - pm, -0.1 - pm, -0.1 - 1i, 0.1 - 1i, 0.1 - pm, -xmax - pm, 999];
            
        case 'mesh_BEDE_up'
            xmax=0.8;
            f = [-xmax + pm, -xmax + 1i, xmax + 1i, xmax + pm, NaN, ...
                xmax - pm, xmax - 1i, -xmax - 1i, -xmax - pm, -0.1 - pm, -0.1 - 1i, 0.1 - 1i, 0.1 - pm, -xmax - pm, NaN, ...
                -xmax + pm, -xmax + 1i, -xmax + pm, 999];
            
            % Extra one
        case 'mesh_BEBE_down'
            xmax=0.1; % Diagonal dodging
            f = [-xmax + pm, -xmax + 1i, xmax + 1i, xmax + pm, NaN, ...
                xmax - pm, xmax - 1i, -xmax - 1i, -xmax - pm, NaN, ...
                -xmax + pm, -xmax + 1i, xmax + 1i, xmax + pm, NaN, ...
                xmax - pm, xmax - 1i, -xmax - 1i, -xmax - pm, -0.9 - pm, 999];
            
        case 'mesh_aAB_down'
            phi0=0.8; psi0=1;
            f = [phi0 + pm,  phi0 + psi0*1i, -phi0 + psi0*1i, -phi0 + pm, NaN, ...
                -phi0 - pm, -phi0 - psi0*1i, phi0 - psi0*1i, phi0 - pm, -phi0 - pm, 999];
            %p = [-phi0 + pm, -phi0 + psi0*1i, phi0 + psi0*1i, phi0 + pm, NaN, ...
            %   phi0 - pm, phi0 - psi0*1i, -phi0 - psi0*1i, -phi0 - pm, 999]
            
        case 'mesh_aABlab_up'
            phi0=0.8; psi0=1;
            f = [-phi0 + pm, -phi0 + psi0*1i, phi0 + psi0*1i, phi0 + pm, NaN, ...
                phi0 - pm, phi0 - psi0*1i, -phi0 - psi0*1i, -phi0 - pm, NaN, ...
                -phi0 + pm, -phi0 + psi0*1i, phi0 + psi0*1i, phi0 + pm, -phi0 + pm, 999];
            
        case 'mesh_aABlab_down'
            phi0=0.8; psi0=1;
            f = [-phi0 + pm, -phi0 + psi0*1i, phi0 + psi0*1i, phi0 + pm, NaN, ...
                phi0 - pm, phi0 - psi0*1i, -phi0 - psi0*1i, -phi0 - pm, NaN, ...
                -phi0 + pm, -phi0 + psi0*1i, phi0 + psi0*1i, phi0 + pm, NaN, ...
                phi0 - pm, -phi0 - pm, 999];
            
        case 'mesh_aAaABlA_down'
            phi0=0.1; psi0=1;
            f = [-phi0 + pm, -phi0 + psi0*1i, phi0 + psi0*1i, phi0 + pm, NaN, ...
                phi0 - pm, phi0 - psi0*1i, -phi0 - psi0*1i, -phi0 - pm, NaN, ...
                -phi0 + pm, -phi0 + psi0*1i, phi0 + psi0*1i, phi0 + pm, NaN, ...
                phi0 - pm, phi0 - 1.5*1i, -0.8 - 1.5*1i, -0.8 - 0.1*1i, 0.1 - 0.1*1i, 0.1 - 0.5*1i, -0.1 - 0.5*1i, -0.1 - pm, 0.5 - pm, 999];

            % ========================================================
            %                  DUAL MESHES
            %
            %               Used for publications
            % * It seems we need at least two non-trivial elements
            %    following a NaN declaration!
            % ========================================================            
        case 'dualmesh.0'
            xmax = 1.2;
            base = [-0.4 + pm, -xmax + pm];
            path.up = [base, 999];
            path.down = [base, NaN, -xmax - pm, -xmax - 0.1i, -xmax - pm, 999];   
            path.f = [-0.4 + pm, -0.4 + 0.1i];
        case 'dualmesh.A'
            xmax = 1.2;
            topmax = 0.3;
            base = [0.4 + pm, 0.4 + topmax*1i, -0.1 + topmax*1i, -0.1 + pm, -xmax + pm];
            path.up = [base, 999];
            path.down = [base, NaN, -xmax - pm, -xmax - 0.05i, -xmax - pm, 999];
            path.f = path.up(1:end-1);
        case 'dualmesh.A.rev'
            xmax = 1.9;
            base = [-0.4 + pm, -0.4 + 1i, 0.1 + 1i, 0.1 + pm, xmax + pm, -xmax + pm];
            path.up = [base, 999];
            path.down = [base, NaN, -xmax - pm, -xmax - 0.05i, -xmax - pm, 999];        
        
        case 'dualmesh.AA_'
            xmax = 1.2;
            base = [0.4 + pm, 0.4 + 1i, -0.1 + 1i, -0.1 + pm, NaN, ...
                    -0.1 - pm, -0.1 - 1i, 0.1 - 1i, 0.1 - pm, ...
                    -xmax - pm];
            path.down = [base, 999];           
            path.up = [base, NaN, -xmax + pm, -xmax + 0.05i, -xmax + pm, 999];
            path.f = path.down(1:end-1);
            
        case 'dualmesh.AA1_A'
            xmax = 1.9;
            base = [0.4 + pm, 0.4 + 1i, -0.4 + 1i, -0.4 + pm, ...
                NaN, -0.4 - pm, -0.4 - 0.1i, 0.9 - 0.1i, ...
                0.9 - 1i, 1.1 - 1i, 1.1 - pm, ...
                NaN, 1.1 + pm, 0.1 + pm, 0.1 + 1i, -0.1 + 1i, ...
                -0.1 + pm, -xmax + pm];
            path.up = [base, 999];
            % For some reason we need to add some fudge elements to end
            path.down = [base, NaN, -xmax - pm, -xmax - 0.05i, -xmax - pm, 999];
            path.f = path.down(1:end-1);
            % ========================================================
            %                  PUBLICATION PATHS
            %
            %               Used for publications
            %
            % ========================================================
        case 'pub_AA_' % for ep = 0.1 (sing about psi = 0.35)
            phi0 = 0.4;
            phi1 = -0.6;
            phi2 = 0.8;
            f = [phi0 + pm, ...
                phi0 + 0.5*1i ...
                phi1 + 0.5*1i, ...
                phi1 + pm, ...
                NaN, ...
                phi1 - pm, ...
                phi1 - 0.65*1i, ...
                phi2 - 0.65*1i, ...
                phi2 - pm, ...
                NaN, ...
                phi2 + pm, ...
                phi2 + 1i];
        case 'pub_topsing_twice' % for ep = 0.1 (sing about psi = 0.35)
            phi0 = 0.4;
            
            f = [phi0 + pm, ...
                phi0 + 0.5*1i ...
                -0.35 + 0.5*1i, ...
                -0.35 + 0.3*1i, ...
                0.3 + 0.3*1i, ...
                0.3 + 0.6*1i ...
                -0.2 + 0.6*1i, ...
                -0.2 + 0.2*1i, ...
                0.4 + 0.2*1i];
        case 'pub_bottomsing_twice'
            phi0 = 0.4;
            f = [phi0 + pm, ...
                phi0 + 0.5*1i ...
                -0.35 + 0.5*1i, ...
                -0.35 + pm, ...
                NaN, ...
                -0.35 - pm, ...
                -0.35 - 0.65*1i, ...
                0.2 - 0.65*1i, ...
                0.2 - 0.2*1i, ...
                -0.1 - 0.2*1i, ...
                -0.1 - 0.5*1i, ...
                0.3 - 0.5*1i, ...
                0.3 - 0.1*1i, ...
                -0.4 - 0.1*1i, ...
                ];
        case 'pub.0.A.ep0.3'
            f_A = 1i*0.175;
            numellipse = 50;
            ellipsefac = 0.2; % For shifting ellipse
            
            f = [0 + pm, 0 + 0.05i];
            
            % Add the spiral
            newp = genspiral(f_A, f(end), numellipse, 2, ellipsefac);
            f = [f, newp];
            
            path.f = f;
            path.singstart = 2;
        case 'pub.A.A_.ep0.3'
            f_A = 1i*0.175;
            numcirc = 50;
            numellipse = 50;
            ellipsefac = 0.2; % For shifting ellipse
            minratuv = 0.8; % For smoothing (around 1)
            nadd = 5; % Number of points to add
            
            
            f = [0 + pm, 0 + 0.05i];
            
            % Add circle
            newp = gencirc(f_A, f(end), numcirc, 1);
            f = [f, newp, 0 + pm];
            
            % Smooth second point and second-to-last point
            f = getSmoothCorner(f, [2, length(f)-1], 'MinRatUV', minratuv, 'NAdd', nadd, 'Display', 'off');
            
            % Cross axis
            f = [f, NaN, 0 - pm, 0 - 0.05i];
            path.singstart = length(f);
            
            % Add the spiral
            newp = genspiral(-f_A, f(end), numellipse, 2, ellipsefac);
            f = [f, newp];
            path.f = f;
        case 'pub.AA_.B'
            f_A = 1i*0.175;
            f_B = 0.788 - 1i*0.554;
            f_B = 0.3 - 1i*0.554;
            
            minratuv = 1; % For smoothing (around 1)
            nadd = 10; % Number of points to add
            numellipse = 50;
            ellipsefac = 0.4; % For shifting ellipse
            
            f = [0.1 + pm, 0.1 + 0.3i, -0.1 + 0.3i];
            f = getSmoothCorner(f, 2, 'MinRatUV', minratuv, 'NAdd', nadd, 'Display', 'off');
            
            f = [f, -0.1 + pm];
            f = getSmoothCorner(f, length(f)-1, 'MinRatUV', minratuv, 'NAdd', nadd, 'Display', 'off');
            
            f = [f, NaN, -0.1 - pm, -0.1 - 0.4i, 0.1 - 0.4i];
            f = getSmoothCorner(f, length(f)-1, 'MinRatUV', minratuv, 'NAdd', nadd, 'Display', 'off');
            
            f = [f, real(f_B) - 0.4i];
            path.singstart = length(f);
            
            % Add the spiral
            newp = genspiral(f_B, f(end), numellipse, 2, ellipsefac);
            f = [f, newp];
            
            path.f = f;
%         case 'pub.AA_.B_ep0.5'
%             f_A = 1i*0.097;
%             f_B = 0.545 - 0.387i;
%             
%             minratuv = 1; % For smoothing (around 1)
%             nadd = 10; % Number of points to add
%             numellipse = 50;
%             ellipsefac = 0.2; % For shifting ellipse
%             
%             phi_r = 0.3;
%             phi_l = -0.3;
%             
%             f = [phi_r + pm, phi_r + 0.3i, phi_l + 0.3i];
%             f = getSmoothCorner(f, 2, 'MinRatUV', minratuv, 'NAdd', nadd, 'Display', 'off');
%             
%             f = [f, phi_l + pm];
%             f = getSmoothCorner(f, length(f)-1, 'MinRatUV', minratuv, 'NAdd', nadd, 'Display', 'off');
%             
%             f = [f, NaN, phi_l - pm, phi_l - 0.25i, 0.1 - 0.25i];
%             f = getSmoothCorner(f, length(f)-1, 'MinRatUV', minratuv, 'NAdd', nadd, 'Display', 'off');
%             
%             f = [f, real(f_B) - 0.25i];
%             path.singstart = length(f);
%             
%             % Add the spiral
%             newp = genspiral(f_B, f(end), numellipse, 2, ellipsefac);
%             f = [f, newp];
%             
%             path.f = f;
        case 'pub.AA_.B_ep0.5'
            f_A = 1i*0.097;
            f_B = 0.545 - 0.387i;
            
            %{
        
            minratuv = 1; % For smoothing (around 1)
            nadd = 10; % Number of points to add
            
            
            phi_r = 0.3;
            phi_l = -0.3;
            
            f = [phi_r + pm, phi_r + 0.3i, phi_l + 0.3i];
            f = getSmoothCorner(f, 2, 'MinRatUV', minratuv, 'NAdd', nadd, 'Display', 'off');
            
            f = [f, phi_l + pm];
            f = getSmoothCorner(f, length(f)-1, 'MinRatUV', minratuv, 'NAdd', nadd, 'Display', 'off');
            
            f = [f, NaN, phi_l - pm, phi_l - 0.25i, 0.1 - 0.25i];
            f = getSmoothCorner(f, length(f)-1, 'MinRatUV', minratuv, 'NAdd', nadd, 'Display', 'off');
            
            f = [f, real(f_B) - 0.25i];
            path.singstart = length(f);
            
            % Add the spiral
            newp = genspiral(f_B, f(end), numellipse, 2, ellipsefac);
            f = [f, newp];
            %}
            
            f = 0.3 + pm;
            newf = genellipse(f(1), 0, 0.387, [0, pi], 30); 
            f = [f, newf];
            
            % Replace last coordinate with jump
            f(end) = real(f(end)) + pm;
            f = [f, NaN, real(f(end)) - pm];
            
            % Note we fudge the next semi-ellipse to use the exact real
            % value so there is no error in the angle
            newf = genellipse(real(f(end)), 0, 0.387, [0, pi/2], 30); 
            f = [f, newf, 0.3 + imag(f_B)*1i];
            
            path.singstart = length(f);
            % Add the spiral
            numellipse = 50;
            ellipsefac = 0.3; % For shifting ellipse
            newp = genspiral(f_B, f(end), numellipse, 2, ellipsefac);
            
            f = [f, newp];
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
                phi0 + pm, phi0 + 1i, -phi0 + 1i, -phi0 + pm, phi0 + pm];
            
        case 'round.A.A_'
            phi0 = 0.1;
            psim = 0.1i;
            f = [phi0 + pm, phi0 + 1i, -phi0 + 1i, -phi0 + pm, ...
                NaN, -phi0 - pm, -phi0 - psim, -phi0 - 1i, phi0 - 1i, phi0 - psim, ...
                -phi0 - psim, -phi0 - 1i, phi0 - 1i, phi0 - psim, ...
                -phi0 - psim];
        case 'round.A.A_.triple' % Triple rotation
            phi0 = 0.1;
            psim = 0.1i;
            f = [phi0 + pm, phi0 + 1i, -phi0 + 1i, -phi0 + pm, ...
                NaN, -phi0 - pm, ...
                -phi0 - psim, -phi0 - 1i, phi0 - 1i, phi0 - psim, ...
                -phi0 - psim, -phi0 - 1i, phi0 - 1i, phi0 - psim, ...
                -phi0 - psim, -phi0 - 1i, phi0 - 1i, phi0 - psim, ...
                -phi0 - psim];
            
            % Need to update notation:
            
        case 'a-0' % Minus sign needed
            phi0 = 0.4;
            
            f = [phi0 + pm, ...
                phi0 + 2*1i ...
                -phi0 + 2*1i, ...
                -phi0 + pm, ...
                phi0 + pm, ...
                phi0 + 2*1i ...
                -phi0 + 2*1i, ...
                -phi0 + pm, ...
                phi0 + pm
                ];
            
            
        case 'A-a' % Minus sign needed
            phi0 = 0.2;
            
            f = [-phi0 + pm, ...
                -phi0 + 1i, ...
                phi0 + 1i, ...
                phi0 + pm, ...
                NaN, ...
                phi0 - pm, ...
                phi0 - 1i, ...
                -phi0 - 1i, ...
                -phi0 - pm, ...
                phi0 - pm, ...
                phi0 - 1i, ...
                -phi0 - 1i, ...
                -phi0 - pm, ...
                phi0 - pm
                ];
            
        case 'B^-aA' % No Minus. Start counting at 8
            C = 1.5; D = 0.1; E = 1.5; F = 0.1;
            f = [-0.4 + pm, -0.4 + 1i, 0.4 + 1i, 0.4 + pm, ...
                NaN, ...
                0.4 - pm, 0.4 - C*1i, -0.1 - C*1i, ...
                -D - C*1i, -E - C*1i, -E - F*1i, -D - F*1i, -D - C*1i , -E - C*1i, -E - F*1i, -D - F*1i, -D - C*1i];
            %{
    case 'B-aA' % Minus
        C = 2; D = 0.1; E = 1.5; F = 0.1;
        p = [-0.4 + pm, -0.4 + 1i, 0.4 + 1i, 0.4 + pm, ...
            NaN, ...
            0.4 - pm, ...
            0.4 - 1i, ...
            -0.4 - 1i, ...
            -0.4 - pm, ...
            0.4 - pm, ...
            0.4 - C*1i, 0.1 - C*1i, ...
            D - C*1i, E - C*1i, E - F*1i, D - F*1i, D - C*1i , E - C*1i, E - F*1i, D - F*1i, D - C*1i];
                %}
                
        case 'B-aA' % Minus
            C=1; D=0.1; E=1; F=0.1;
            f = [-0.4 + pm, -0.4 + 1i, 0.4 + 1i, 0.4 + pm, NaN, ...
                0.4 - pm, 0.4 - 1i, -D - 1i, ...
                -D - F*1i, -C - F*1i, -C - E*1i, -D - E*1i, -D - F*1i, -D - F*1i, -C - F*1i, -C - E*1i, -D - E*1i, -D+F*1i];
            
            
        case 'C^-aAaA' % No Minus
            phi0 = 0.2;
            psi0 = 1;
            C = 1.5; D = 0.1; E = 2.5; F = 0.1;
            f = [-phi0 + pm, -phi0 + psi0*1i, phi0 + psi0*1i, phi0 + pm, ...
                NaN, ...
                phi0 - pm, phi0 - psi0*1i, -phi0 - psi0*1i, -phi0 - pm, ...
                NaN, ...
                -phi0 + pm, -phi0 + psi0*1i, phi0 + psi0*1i, phi0 + pm, ...
                NaN, ...
                phi0 - pm, ...
                phi0 - C*1i, ...
                -D - C*1i, -E - C*1i, -E - F*1i, -D - F*1i, -D - C*1i , -E - C*1i, -E - F*1i, -D - F*1i, -D - C*1i];
            
        case 'C-aAaA' % Minus
            phi0 = 0.2;
            psi0 = 1;
            C = 1.5; D = 0.1; E = 1.5; F = 0.1;
            f = [-phi0 + pm, -phi0 + psi0*1i, phi0 + psi0*1i, phi0 + pm, ...
                NaN, ...
                phi0 - pm, phi0 - psi0*1i, -phi0 - psi0*1i, -phi0 - pm, ...
                NaN, ...
                -phi0 + pm, -phi0 + psi0*1i, phi0 + psi0*1i, phi0 + pm, ...
                NaN, ...
                phi0 - pm, ...
                phi0 - C*1i, ...
                -D - C*1i, -E - C*1i, -E - F*1i, -D - F*1i,  -D - C*1i, -E - C*1i, -E - F*1i, -D - F*1i, -D - C*1i];
            
        case 'a-aA'
            phi0 = 0.1; psi0 = 1;
            f = [-phi0 + pm, -phi0 + psi0*1i, phi0 + psi0*1i, phi0 + pm, NaN, ...
                phi0 - pm, phi0 - psi0*1i, -phi0 - psi0*1i, -phi0 - pm, NaN, ...
                -phi0 + pm, -phi0 + psi0*1i, phi0 + psi0*1i, phi0 + pm, -phi0 + pm, -phi0 + psi0*1i, phi0 + psi0*1i, phi0 + pm, -phi0 + pm];
            
        case 'b-aAB^a' %13
            phi0 = 0.1; psi0 = 1;
            f = [-phi0 + pm, -phi0 + psi0*1i, phi0 + psi0*1i, phi0 + pm, NaN, ...
                phi0 - pm, phi0 - psi0*1i, -0.8 - psi0*1i, -0.8 - pm, NaN, ...
                -0.8 + pm, -0.8 + psi0*1i, ...
                0.1 + psi0*1i, 1.5 + psi0*1i, 1.5 + pm, 0.1 + pm, 0.1 + psi0*1i, 1.5 + psi0*1i, 1.5 + pm, 0.1 + pm, 0.1 + psi0*1i];
            
        case 'a-aAB^' %12
            phi0 = 0.7; psi0 = 1;
            f = [-phi0 + pm, -phi0 + psi0*1i, phi0 + psi0*1i, phi0 + pm, NaN, ...
                phi0 - pm, phi0 - psi0*1i, -phi0 - psi0*1i, -phi0 - pm, NaN, ...
                -phi0 + pm, -0.1 + pm, 0.1 + pm, 0.1 + 0.3*1i, -0.1 + 0.3*1i, -0.1 + pm, 0.1 + pm, 0.1 + 0.3*1i -0.1 + 0.3*1i, -0.1 + pm];
            
        case 'B-aAB^' %13
            phi0 = 0.9; psi0 = 1; A = 0.8; B=0.8; C=0.1; D=0.01;
            f = [-phi0 + pm, -phi0 + psi0*1i, phi0 + psi0*1i, phi0 + pm, NaN, ...
                phi0 - pm, phi0 - psi0*1i, -phi0 - psi0*1i, -phi0 - pm, NaN, ...
                -phi0 + pm, -phi0 + D*1i, A + D*1i, A + B*1i, C + B*1i, C + D*1i, A + D*1i,  A + B*1i, C + B*1i, C + D*1i, A + D*1i];
            
            
            
            
        case 'left_diag_singularity3'
            phi0 = 0.2;
            C = 3; D = 0.1; E = 2.5; F = 0.05;
            f = [-phi0 + pm, -phi0 + 1i, phi0 + 1i, phi0 + pm, ...
                NaN, ...
                phi0 - pm, ...
                phi0 - 1i, ...
                -phi0 - 1i, ...
                -phi0 - pm, ...
                NaN, ...
                -phi0 + pm, ...
                -phi0 + 1i, ...
                phi0 + 1i, ...
                phi0 + pm, ...
                NaN, ...
                phi0 - pm, ...
                phi0 - 1i, ...
                -phi0 - 1i, ...
                -phi0 - pm, ...
                NaN, ...
                -phi0 + pm, ...
                -phi0 + 1i, ...
                phi0 + 1i, ...
                phi0 + pm, ...
                NaN, ...
                phi0 - pm, ...
                phi0 - C*1i, -0.1 - C*1i, ...
                -D - C*1i, -E - C*1i, -E - F*1i, -D - F*1i, -D - C*1i , -E - C*1i, -E - F*1i, -D - F*1i, -D - C*1i];
            
            % ========================================================
            %                   CIRCLE PATHS
            %
            %    Same as round paths but with circles
            %
            % ========================================================
        case 'circle.0.A.ep0.3'
            numtours = 2;
            numcirc = 20;
            z0 = 1i*0.1760459658;
            
            f = [0 + pm, 0 + 0.05i];
            newp = gencirc(z0, f(end), numcirc, numtours);
            f = [f, newp];
            path.f = f;
            path.singstart = 3;
        case 'circle.A.A_.ep0.3'
            f_A = 1i*0.175;
            %             numcirc = 20;
            numcirc = 50; % Pub quality?
            
            f = [0 + pm, 0 + 0.05i];
            newp = gencirc(f_A, f(end), numcirc, 1);
            f = [f, newp, 0 + pm, NaN, 0 - pm, 0 - 0.05i];
            newp = genspiral(-f_A, f(end), numcirc, 2, 0.2);
            
            path.singstart = length(f) + 1;
            
            path.f = [f, newp];            
            
        case 'circle.AA_.B.ep0.3' % Centred on 0.3 but reasonable for ep \in [0.2, 0.5]
            f_B = 0.788 - 1i*0.554;
            numcirc = 10;
            
            f = [0.1 + pm, 0.1 + 0.3*1i, -0.1 + 0.3*1i, -0.1 + pm, NaN, ...
                -0.1 - pm, -0.1 - 1i*0.5, 0.3 - 1i*0.554];
            
            newp = gencirc(f_B, f(end), numcirc, 2);
            
            path.singstart = 8;
            path.f = [f, newp];
            
        case 'circle.AA_.B.ep0.2' % Centred on 0.2 but reasonable for ep \in ?
            f_B = 0.97 - 1i*0.65;
            numcirc = 10;
            
            f = [0.1 + pm, 0.1 + 0.5*1i, -0.1 + 0.5*1i, -0.1 + pm, NaN, ...
                -0.1 - pm, -0.1 - 0.65*1i, 0.6 - 0.65*1i];
            
            newp = gencirc(f_B, f(end), numcirc, 2);
            
            path.singstart = 8;
            path.f = [f, newp];
            
        case 'circle.AA_.B.ep0.1' % Centred on 0.1 but reasonable for ep \in ?
            f_B = 1.26 - 1i*0.8;
            numcirc = 10;
            
            f = [0.1 + pm, 0.1 + 0.5*1i, -0.1 + 0.5*1i, -0.1 + pm, NaN, ...
                -0.1 - pm, -0.1 - 0.8*1i, 0.85 - 0.8*1i];
            
            newp = gencirc(f_B, f(end), numcirc, 2);
            
            path.singstart = 8;
            path.f = [f, newp];
           
        case 'circle.AA_.Bl.ep0.05' % Centred on 0.05 but reasonable for ep \in ?
            f_B = -1.54 - 1i*0.94;
            numcirc = 20;
            
            f = [-0.1 + pm, -0.1 + 0.8*1i, 0.1 + 0.8*1i, 0.1 + pm, NaN, ...
                0.1 - pm, 0.1 - 1i*0.94, -0.8 - 1i*0.94];
            
            newp = gencirc(f_B, f(end), numcirc, 2);
            
            path.singstart = 8;
            path.f = [f, newp];

            
        case 'circle.AA_.Bl.ep0.1' % Centred on 0.1 but reasonable for ep \in [0.2, 0.5]
            f_B = -1.26 - 1i*0.8;
            numcirc = 20;
            
            f = [-0.1 + pm, -0.1 + 0.5*1i, 0.1 + 0.5*1i, 0.1 + pm, NaN, ...
                0.1 - pm, 0.1 - 1i*0.8, -0.9 - 1i*0.8];
            
            newp = gencirc(f_B, f(end), numcirc, 2);
            
            path.singstart = 8;
            path.f = [f, newp];
            
        case 'circle.AA_.Bl.ep0.3' % Centred on 0.3 but reasonable for ep \in [0.2, 0.5]
            f_B = -0.788 - 1i*0.554;
            numcirc = 20;
            phimax = 0.5;
            psimax = 0.5;
            
            
            f = [-phimax + pm, -phimax + psimax*1i, phimax + psimax*1i, phimax + pm, NaN, ...
                phimax - pm, phimax - 1i*0.554, -psimax - 1i*0.554];
            
            newp = gencirc(f_B, f(end), numcirc, 2);
            
            path.singstart = 8;
            path.f = [f, newp];
            
        case 'circle.AA_.Bl.ep0.5' % Centred on 0.5 but reasonable for ep \in ?
            f_B = -0.53 - 1i*0.40;
            numcirc = 20;
            
            f = [-0.1 + pm, -0.1 + 0.3*1i, 0.1 + 0.3*1i, 0.1 + pm, NaN, ...
                0.1 - pm, 0.1 - 1i*0.4, -0.4 - 1i*0.4];
            
            newp = gencirc(f_B, f(end), numcirc, 2);
            
            path.singstart = 8;
            path.f = [f, newp];
        
        case 'circle.AA_.Bl.ep0.6' % Centred on 0.6 but reasonable for ep \in ?
            f_B = -0.43 - 1i*0.33;
            numcirc = 20;
            
            f = [-0.1 + pm, -0.1 + 0.3*1i, 0.1 + 0.3*1i, 0.1 + pm, NaN, ...
                0.1 - pm, 0.1 - 1i*0.33, -0.3 - 1i*0.33];
            
            newp = gencirc(f_B, f(end), numcirc, 2);
            
            path.singstart = 8;
            path.f = [f, newp];
            
        case 'circle.AA_.Bl.ep0.7' % Centred on 0.7 but reasonable for ep \in ?
            f_B = -0.34 - 1i*0.25;
            numcirc = 20;
            
            f = [-0.1 + pm, -0.1 + 0.3*1i, 0.1 + 0.3*1i, 0.1 + pm, NaN, ...
                0.1 - pm, 0.1 - 1i*0.25, -0.2 - 1i*0.25];
            
            newp = gencirc(f_B, f(end), numcirc, 2);
            
            path.singstart = 8;
            path.f = [f, newp];
            
        case 'circle.AA_.Bl.ep0.8' % Centred on 0.8 but reasonable for ep \in ?
            f_B = -0.23 - 1i*0.16;
            numcirc = 20;
            
            f = [-0.1 + pm, -0.1 + 0.3*1i, 0.1 + 0.3*1i, 0.1 + pm, NaN, ...
                0.1 - pm, 0.1 - 1i*0.16, -0.175 - 1i*0.16];
            
            newp = gencirc(f_B, f(end), numcirc, 2);
            
            path.singstart = 8;
            path.f = [f, newp];    
            
            
        case 'circle.AA_AA_.Cl.ep0.4' % Centred on 0.4 but reasonable for ep \in ?
            f_C = -0.98 - 0.33*1i;
            numcirc = 20;
            
            f = [-0.1 + pm, -0.1 + 0.3*1i, 0.1 + 0.3*1i, 0.1 + pm, NaN, ...
                0.1 - pm, 0.1 - 0.3*1i, -0.1 - 0.3*1i, -0.1 - pm, NaN, ...
                -0.1 + pm, -0.1 + 0.3*1i, 0.1 + 0.3*1i, 0.1 + pm, NaN, ...
                0.1 - pm, 0.1 - 0.6*1i, -0.7 - 1i*0.6, -0.7 - 0.33*1i];
            
            newp = gencirc(f_C, f(end), numcirc, 2);
            
            path.singstart = 19;
            path.f = [f, newp];
            
        case 'circle.AA_AA_.C.ep0.4' % Centred on 0.4 but reasonable for ep \in ?
            f_C = 0.91 - 0.41*1i;
            numcirc = 20;
            psimax = 0.5;
            
            
            f = [0.1 + pm, 0.1 + psimax*1i, -0.1 + psimax*1i, -0.1 + pm, NaN, ...
                -0.1 - pm, -0.1 - psimax*1i, 0.1 - psimax*1i, 0.1 - pm, NaN, ...
                0.1 + pm, 0.1 + psimax*1i, -0.1 + psimax*1i, -0.1 + pm, NaN, ...
                -0.1 - psimax*1i, 0.6 - psimax*1i, 0.6 - 0.43*1i];
            
            newp = gencirc(f_C, f(end), numcirc, 2);
            
            path.singstart = length(f);
            path.f = [f, newp];
            
        case 'circle.AA_AA_.C.ep0.2' % Centred on 0.2 but reasonable for ep \in ?
            %f_C = 0.91 - 0.69*1i;
            %f_C = 0.82 - 0.97*1i;
            f_C = 1.21 - 1.07*1i;
            numcirc = 20;
            psimax = 0.8;
            phimax = 0.8;
            
            f = [phimax + pm, phimax + psimax*1i, -phimax + psimax*1i, -phimax + pm, NaN, ...
                -phimax - pm, -phimax - psimax*1i, phimax - psimax*1i, phimax - pm, NaN, ...
                phimax + pm, phimax + psimax*1i, -phimax + psimax*1i, -phimax + pm, NaN, ...
                -phimax - psimax*1i, 0.7 - psimax*1i, 0.7 - 1.07*1i];
            
            
            newp = gencirc(f_C, f(end), numcirc, 2);
            
            path.singstart = length(f);
            path.f = [f, newp];  
           
            
        case 'circle.AA_AA_.C.ep0.1' % Centred on 0.1 but reasonable for ep \in ?
            f_C = 0.90 - 0.79*1i;
            numcirc = 20;
            psimax = 1;
            phimax = 0.05;
            
            f = [phimax + pm, phimax + psimax*1i, -phimax + psimax*1i, -phimax + pm, NaN, ...
                -phimax - pm, -phimax - psimax*1i, phimax - psimax*1i, phimax - pm, NaN, ...
                phimax + pm, phimax + psimax*1i, -phimax + psimax*1i, -phimax + pm, NaN, ...
                -phimax - psimax*1i, 0.7 - psimax*1i, 0.7 - 0.79*1i];
            
            newp = gencirc(f_C, f(end), numcirc, 2);
            
            path.singstart = 18;
            path.f = [f, newp];                
                                  
                        
            
        otherwise
            warning('Could not find path')
            keyboard
    end
end % end main function

function z = gencirc(z0, z1, numpts, numtours)
    % GENCIRC generates points around a circle
    rad = abs(z1 - z0);
    theta0 = angle(z1 - z0);
    t = linspace(0, numtours*2*pi, numtours*numpts + 1);
    % Remove the first point, which repeats
    t = t(2:end);
    z = z0 + rad.*exp(1i*(theta0 + t));
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