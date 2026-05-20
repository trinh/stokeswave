function patchdat = getPatchForComplexContour(C, ratio_normal, reimname)
% DRAW_PATCHES draws a series of patches along a curve

projectxy = 'true';

plotme = 0;
if nargin == 0
    plotme = 1;
    
    ratio_normal = 1;
    
    s = linspace(0, 2*pi, 50);
    a = 0.5;
        
    gam = (s+a).*exp(1i*s);
    f = gam.^(1/2);
    
    % Use 'real', 'imag', 'abs', etc.
    reimname = 'real';
    C = [gam(:).'; f(:).'];
end

reim = str2func(reimname);

if plotme == 1
    xx = linspace(-10, 10, 50);
    yy = xx;
    [XX,YY] = meshgrid(xx,yy);
    ft = @(z) z.^(1/2);
    ZZ = real(ft(XX + 1i*YY));
    figure
    surf(XX, YY,ZZ);
    hold on
end


%Use a 2nd order approximation for the derivatives of the trajectory
if size(C,2) >= 3 
    dC = [C(:,1:3)*[-3; 4; -1]/2  [C(:,3:end) - C(:,1:end-2)]/2 C(:,end-2:end)*[1; -4; 3]/2];   
else
    dC = C(:,[2 2]) - C(:,[1 1]);
end
% Because of the conformality, this is the 90 degree normal vector
idC = -1i*dC;

N_contour = size(dC,2);


verts_right = [];
faces_right = [];
verts_left = [];
faces_left = [];
for j = 1:N_contour - 1
    
    % if isnan(gam(j)) == 1 || isnan(gam(j+1)) == 1
    %    continue;
    % end
    
    % Find the normals to the 2D curve, and then use the derivative of the
    % surface function to calculate height raise
    n0 = [real(idC(1,j)), imag(idC(1,j)), reim(idC(2,j))];
    n0 = n0/norm(n0);    
    
    n1 = [real(idC(1,j+1)), imag(idC(1,j+1)), reim(idC(2,j+1))];
    n1 = n1/norm(n1);
    
    % Actually I want to scale the normals so that the projected distance
    % in the x-y plane is ratio_normal
    switch projectxy
        case 'true'
            n0_lengthxy = norm(n0(1:2));
            n1_lengthxy = norm(n1(1:2));
            n0 = n0/n0_lengthxy;
            n1 = n1/n1_lengthxy;
    end

    A1 = [real(C(1,j)), imag(C(1,j)), reim(C(2,j))];
    A4 = [real(C(1,j+1)), imag(C(1,j+1)), reim(C(2,j+1))];
    
    % 4 vertices for each face and N - 1 faces
    verts_right(end+1, :) = A1; 
    verts_right(end+1, :) = A1 + ratio_normal*n0;   % A2
    verts_right(end+1, :) = A4 + ratio_normal*n1;   % A3
    verts_right(end+1, :) = A4;
    verts_left(end+1, :) = A1; 
    verts_left(end+1, :) = A1 - ratio_normal*n0;    % A2_left
    verts_left(end+1, :) = A4 - ratio_normal*n1;    % A3_left
    verts_left(end+1, :) = A4;
    
    % Face index just goes [1, 2, 3, 4; 5 6 7 8; ...]
    faces_right(end+1,:) = 4*(j-1) + [1, 2, 3, 4];
    faces_left(end+1,:) = 4*(j-1) + [1, 2, 3, 4];
    
end
if nargin == 0
    p = patch('Vertices', verts_right, 'Faces', faces_right);
    set(p, 'FaceColor', 'w');
    hold off
end

patchdat = [];
patchdat.verts_right = verts_right;
patchdat.faces_right = faces_right;
patchdat.verts_left = verts_left;
patchdat.faces_left = faces_left;