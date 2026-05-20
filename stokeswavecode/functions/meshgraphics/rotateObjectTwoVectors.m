function Cnew = rotateObjectTwoVectors(C, u, v)
% rotateObjectTwoVectors rotates an object so that vector u points towards
% vector v
%   C = 3 x N matrix


u = u(:)/norm(u);
v = v(:)/norm(v);

% The axis of rotation is perpendicular u and v
z = cross(u, v);

if norm(z) ~= 0
    % Original
    z = z/norm(z);
    
    % Calculate the angle between old tangent and new tangent
    theta = real(acos(dot(u, v)));
    
    Cnew = rotateObjectAngleAxis(C, theta, z);
end



% Below is more efficient but not as clean
%{
% We want to rotate [0;0;1] 180degrees around an axis 'z' to become dC
u = u(:)/norm(u);
v = v(:)/norm(v);

npt = size(C,2);

% The axis of rotation is perpendicular u and v
z = cross(u, v);

if norm(z) ~= 0
    % Original
    z = z/norm(z);
    
    % Calculate the angle between old tangent and new tangent
    q = real(acos(dot(u, v)));
    
    Z = repmat(z,1,npt);
    
    % Rodrigues' formula    
    Cnew = C*cos(q) + cross(Z,C)*sin(q)+Z*(1-cos(q))*diag(dot(Z,C));
%     Cnew2 = rotateObjectAngleAxis(C, q, z);
%     keyboard
end
%}
