clear

close all

path = [-1 - 1i, ...
        1 - 1i, ...
        1 + 1i ...
        -1 + 1i ...
        -1 - 1i];
    
z0 = 1 + 1i;
rad = 1;
path = [z0 - rad - rad*1i, ...
        z0 + rad - rad*1i, ...
        z0 + rad + rad*1i, ...
        z0 - rad - rad*1i, ...
        z0 - rad - rad*1i];
    
    
[Z, Y] = odecontour(path);


% % Scheme to mesh from the t-plane
% theta = linspace(0, 2*pi, 20);
% r = linspace(1, 3, 100);
% [R, T] = meshgrid(r, theta);
% 
% solmat = 0*R;
% for j = 1:length(theta)
%     solmat(j, :) = ode45(fwd, r, initcond);
% end
% % Convert from Polar to Cartesian
% [X, Y] = pol2cart(T, R);
% % Plot surf(X, Y, solmat)