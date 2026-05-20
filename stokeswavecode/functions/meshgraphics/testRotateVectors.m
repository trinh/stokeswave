% Testing rotations using the rotateObjectTwoVectors code

clear
close all

% Create an object
[x,y,z] = peaks(20);
% We wish to rotate such that u is sent to v
u = [0, 0, 1];
v = [0, 1, 0];
for j = 1:size(x,1)
    for k = 1:size(x,2)
        tmp = rotateObjectTwoVectors([x(j,k); y(j,k); z(j,k)], u, v);
        xnew(j,k) = tmp(1);
        ynew(j,k) = tmp(2);
        znew(j,k) = tmp(3);
    end
end

surf(x, y, z, 'FaceColor', 'b', 'FaceAlpha', 0.5); hold on;
surf(xnew, ynew, znew, 'FaceColor', 'r', 'FaceAlpha', 0.5)
xlabel('x');
ylabel('y');
zlabel('z');