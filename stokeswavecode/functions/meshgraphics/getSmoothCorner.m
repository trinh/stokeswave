% function newpath = getSmoothCorner(path, ind, minr, n_add)
function newpath = getSmoothCorner(varargin)
% getSmoothCorner smooths out a complex contour
%
%   path = (1 x N) complex contour
%   ind = set of locations where you want corners smoothed
%   minr = radius of curvature (try 0.2 for line segments 1)
%   n_add = number of points to add to smooth


% Parse the inputs
myparse = inputParser;
myparse.addOptional('path', [0, 1i, 2], @(x)length(x) > 0);
myparse.addOptional('ind', 2, @(x) length(x) > 0);
% myparse.addParameter('MinR', Inf, @isnumeric);
myparse.addParameter('MinRatUV', 1, @isnumeric);
myparse.addParameter('NAdd', 10, @isnumeric);
myparse.addParameter('Display', 'off', @ischar);
myparse.parse(varargin{:})

path = myparse.Results.path;
ind = myparse.Results.ind;
% minr = myparse.Results.MinR;
n_add = myparse.Results.NAdd;
minratuv = myparse.Results.MinRatUV;

newpath = [];
for j = 1:length(ind)
    minipath = [path(ind(j) - 1), path(ind(j)), path(ind(j)+1)];
    C = [real(minipath(:).'); imag(minipath(:).')];
    
    u = C(:,2) - C(:,1);
    v = C(:,3) - C(:,2);
        
    % Determine intersections with a circle of radius r
    % Take either minr or a proportion of the segment length
%     r = min([minr, norm(u)*minratuv, norm(v)*minratuv]);    
%     r = minr;
    
%     p1 = C(:,2) - r*u/norm(u);
%     p2 = C(:,2) + r*v/norm(v);
%     
%     
%     q1 = C(:,2);
%     q2 = C(:,2);
    p1 = C(:,1);
    q1 = C(:,1) + minratuv*u;
    q2 = C(:,3) - minratuv*v;
    p2 = C(:,3);
    
    % Q = bezierInterp(p1, q1, q2, p2);
    pts = [p1, q1, q2, p2];
    Q = getBezier4(pts, n_add);
    % Note that Q includes the flanking points
    Q = Q(:,2:end-1);
    
    % Replace the path(ind(j)) point with Q
    if j == 1
        newpath = [path(1:ind(j)-1), Q(1,:) + 1i*Q(2,:)];
    else
        newpath = [newpath, path(ind(j-1) + 1:ind(j) - 1), Q(1,:) + 1i*Q(2,:)];
    end

%     plot(real(minipath), imag(minipath)); hold on
%     plot(p1(1), p1(2), 'ro');
%     plot(p2(1), p2(2), 'bo');
%     keyboard
end
newpath = [newpath, path(ind(end) + 1:end)];

switch myparse.Results.Display
    case 'on'
        figure(1) 
        plot(real(path), imag(path), 'b'); hold on
        plot(real(newpath), imag(newpath), 'r.-');
        hold off
end
