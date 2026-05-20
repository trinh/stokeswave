function Cnew = rotateObjectAngleAxis(C, theta, u)

if nargin == 0
    close all
    
    x = linspace(-1, 1, 10);
    y = linspace(-2, 2, 20);
    [X,Y] = meshgrid(x, y);
    Z = sin(X+Y);
    [Ny, Nx] = size(X);
    
    u = [0, 1, 0];
    theta = pi/2;
    
    C = [X(:).'; Y(:).'; Z(:).'];
end

u = u/norm(u);
R = getRotationMatrixAngleAxis(theta,u);

Cnew = R*C;

if nargin == 0
    Xnew = reshape(Cnew(1,:), Ny, Nx);
    Ynew = reshape(Cnew(2,:), Ny, Nx);
    Znew = reshape(Cnew(3,:), Ny, Nx);
       
    surf(X, Y, Z, 'FaceColor', 'b', 'FaceAlpha', 0.5); hold on;
    surf(Xnew, Ynew, Znew, 'FaceColor', 'r', 'FaceAlpha', 0.5)
    xlabel('x');
    ylabel('y');
    zlabel('z');
end