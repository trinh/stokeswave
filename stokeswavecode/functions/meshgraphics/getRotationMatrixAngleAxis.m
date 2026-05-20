function R = getRotationMatrixAngleAxis(theta,u)
% getRotationMatrixAngleAxis The Rodrigues' formula for rotation matrices.
% 
%       ARGUMENT DESCRIPTION:
%           THETA - angle of rotation (radians).
%               U - unit vector
% 
%       OUTPUT DESCRIPTION:
%               R - rotation matrix.
% 

u = u./norm(u,2);
S = [    0  -u(3) u(2);
      u(3)   0   -u(1);
      -u(2) u(1)   0  ];
R = eye(3) + sin(theta)*S + (1-cos(theta))*S^2;