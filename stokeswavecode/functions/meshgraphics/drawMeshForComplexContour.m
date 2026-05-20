function [X, Y, Z] = drawMeshForComplexContour(base, C, reimname)
% drawMeshForComplexContour creates a mesh along a trajectory 
%
%   base = [x;y]    A 2-row Matrix defining the 2d base
%   C = [x + iy; f] A 2-row Matrix defining the 3d path (3xn matrix)
%   reimname        A string 'real', 'imag', 'abs', etc.

alignnormal = 1; % Turn on alignment with normals
plotme = 0;

if nargin < 3
    plotme = 1; % Turn on to plot
    % base = [-1 -1 1 1; 0 -.2 -.2 0]/10;
    base = [-1 -1 1 1; .2 0 0 0.2]/4;

    % Create a test curve
    s = linspace(0, 2*pi, 101);
    a = 0.5;
    gam = 0.5*(s+a).*exp(1i*s);
    
    % Piecewise linear path
    % gam = -2 + 1i*linspace(0, 2, 50);
    % gam = [gam, -2 + 2i + linspace(0, 2, 50)];

    ft = @(z) sqrt(z);
    f = ft(gam);

    C = [gam(:).'; f(:).'];
    % C = 1*C./max(abs(C(:))); %normalize it
    
    % Use 'real', 'imag', 'abs', etc.
    reimname = 'real';
end


reim = str2func(reimname);

if plotme == 1
    xx = linspace(-4, 4, 50); yy = xx;
    [XX,YY] = meshgrid(xx,yy);

    ZZ = reim(ft(XX + 1i*YY));
    figure
    surf(XX, YY, ZZ-0.01, 'FaceColor', 'b', 'FaceAlpha', 0.5); hold on;
end

%% Calculate derivatives and allocate matrices
npt = size(base,2);
base = [base; zeros(1,npt)];

%Use a 2nd order approximation for the derivatives of the trajectory
if size(C,2) >= 3 
    dC = [C(:,1:3)*[-3; 4; -1]/2  [C(:,3:end) - C(:,1:end-2)]/2 C(:,end-2:end)*[1; -4; 3]/2];   
else
    dC = C(:,[2 2]) - C(:,[1 1]);
end

dC0 = find(sum(abs(dC),1) == 0,1);    %Check for stagnation points
if ~isempty(dC0)
    warning('Removing stagnation points found in trajectory');
    dCgood = find(sum(abs(dC),1) ~= 0);
    C = C(:,dCgood); %
    dC = dC(:,dCgood); %
end

N_contour = size(dC,2);
SUR = nan(3,npt,N_contour);

%% Generate and plot the surface

% We define an initial forward direction
dCvec_prev = [0; 0; 1];
% And an initial normal direction
norm_prev = [0; 1; 0];

% Because of the conformality, this is the 90 degree normal vector
idC = -1i*dC;

% We calculate the normals along the contour by crossing the tangential
% vector with the normal vector
v1 = [real(idC(1,:)); imag(idC(1,:)); reim(idC(2,:))];
v2 = [real(dC(1,:)); imag(dC(1,:)); reim(dC(2,:))];
nv = cross(v1, v2);
for k = 1:N_contour
    nv(:,k) = nv(:,k)/norm(nv(:,k));
end



if plotme == 1
    plot3(real(C(1,:)), imag(C(1,:)), reim(C(2,:)), 'k'); hold all
    plot3(base(1,:), base(2,:), base(3,:), 'o-'); 
end
    
% We run through each step in the contour
for k = 1:N_contour
    
    % Define a next forward vector
    dCvec = [real(dC(1,k)); imag(dC(1,k)); reim(dC(2,k))];
    dCvec = dCvec/norm(dCvec);
    % Rotate the current forward vector to the next
    base = rotateObjectTwoVectors(base, dCvec_prev, dCvec);    
    % Rotate also the current normal 
    norm_prev = rotateObjectTwoVectors(norm_prev, dCvec_prev, dCvec); 
    % Set the new prevous forward vector
    dCvec_prev = dCvec;
    if alignnormal == 1
        % Rotate the old normal to the next
        base = rotateObjectTwoVectors(base, norm_prev, nv(:,k));
        % Now we set the current forward vector; note require further
        % rotation due to the change in normals
        dCvec_prev = rotateObjectTwoVectors(dCvec_prev, norm_prev, nv(:,k));
    end
    norm_prev = nv(:,k);
    
    % Should be perpendicular
    % dot(norm_prev, dCvec_prev)
    
    
    SUR(:,:,k) = base + repmat([real(C(1,k)); imag(C(1,k)); reim(C(2,k))], 1, npt);
    if plotme == 1
        pt = [real(C(1,k)); imag(C(1,k)); reim(C(2,k))];
        tmn = [pt, pt + norm_prev];
        plot3(tmn(1,:), tmn(2,:), tmn(3,:), 'r', 'LineWidth', 1);
        test_tang = [pt, pt + dCvec];
        plot3(test_tang(1,:), test_tang(2,:), test_tang(3,:), 'g', 'LineWidth', 1);
    end
end


X = squeeze(SUR(1,:,:));
Y = squeeze(SUR(2,:,:));
Z = squeeze(SUR(3,:,:));

% Fix any NaN's (This should not happen...)
xisnan = find(isnan(X(1,:)), 1);
xnotnan = find(~isnan(X(1,:)));
if ~isempty(xisnan)
    warning('NaN''s found');
    X = X(:,xnotnan);
    Y = Y(:,xnotnan);
    Z = Z(:,xnotnan);
    C = C(:,xnotnan);
end

if plotme == 1
    surf(X,Y,Z)
    axis equal
    dragzoom;
end



