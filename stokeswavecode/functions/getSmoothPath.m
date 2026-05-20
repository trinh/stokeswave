function newpath = getSmoothPath(path, rmin, n_smooth)
% getSmoothPath smooths out a path (a 1xn complex-valued vector)
%   It assumes that all interior points are corners, except those marking a
%   NaN (crossing the axis), for which the flanking points should be
%   removed from the smoothing as well

% Remove any NaN by finding locations
ind_nan = find(isnan(abs(path))==1);
if isempty(ind_nan)
    ind_smooth = 2:length(path)-1;
else
    for j = 1:length(ind_nan)
        if j == 1
            ind_smooth = (1:ind_nan(1) - 2);
        else
            ind_smooth = [ind_smooth, ind_nan(j-1)+2, ind_nan(j) - 1];
        end
    end
    ind_smooth = [ind_smooth, ind_nan(j) + 2:length(path)];
    ind_smooth = ind_smooth(2:end-1);
end
newpath = getSmoothCorner(path, ind_smooth, rmin, n_smooth);    
