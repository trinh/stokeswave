function handles = qt_bdint_path_plotter(path, pathdat, varargin)
% QT_BDINT_PATH_PLOTTER plots a mesh along a path 

% Parse the inputs
myparse = inputParser;
myparse.addRequired('path', @(x) length(x) > 0);
myparse.addRequired('pathdat', @(x) true);
myparse.addParameter('ZPlot', 'real', @ischar)
myparse.addParameter('ZFunc', 'omega', @ischar)
myparse.addParameter('NumArc', 300, @isnumeric)
myparse.addParameter('PatchWidth', 0.01, @isnumeric)
myparse.parse(path, pathdat, varargin{:})

reimname = myparse.Results.ZPlot;
NumArc = myparse.Results.NumArc;
PatchWidth = myparse.Results.PatchWidth;

switch myparse.Results.ZFunc
    case 'omega'
        C = [pathdat.f(:).'; pathdat.omega(:).'];
    case 'dfdz'
        C = [pathdat.f(:).'; exp(pathdat.omega(:).')];
    case 'z'
        C = [pathdat.f(:).'; exp(pathdat.z(:).')];
    otherwise
        warning('wrong ZFunc')
        keyboard
end
Cnew = getComplexCurveByArcLength(C, NumArc, reimname);
handles.surf = drawPatchForComplexContour(Cnew, PatchWidth, reimname);

handles.arcplot = plot_arcprofile(pathdat.f, pathdat.z, path);

%{
%% Universal constants
smoothing = 'on';
reimname = 'real';
reim = str2func(reimname);
figsize = [600 400 600 1000];
myeps = 1e-2;
minz = -10;


%% Parameters
% =====================================
%            Path input
path = get_path('pub_topsing_twice');
rmin = 0.2;
n_smooth = 3;
n_arc = 200;
patch_width = 0.04;


%%
% Smooth the path
switch smoothing 
    case 'on'
        path = getSmoothPath(path, rmin, n_smooth);
end

% Solve for the path
[pathdat,singloc] = qt_bdint_path_solve(path, solseries, []);

% Plot
omega = pathdat.tau - 1i*pathdat.theta;
W = exp(omega);


figure(5); figshift; 
plot3(real(pathdat.f),imag(pathdat.f), reim(W)+myeps, '-', 'LineWidth', 2); hold all
plot3(real(pathdat.f),imag(pathdat.f), 0*reim(W) + minz, '-', 'LineWidth', 2);
xlabel('\phi'); ylabel('\psi');

C = [pathdat.f(:).'; W(:).'];
Cnew = getComplexCurveByArcLength(C, n_arc, reimname);
C = Cnew;
handles = drawPatchForComplexContour(C, patch_width, reimname);

set(gcf, 'Units', 'pixels', 'Position', figsize);
view([-50, 40]);
shg

%%

% =====================================
%            Path input
path = get_path('pub_bottomsing_twice');

%% Parameters
rmin = 0.2;
n_smooth = 3;
n_arc = 200;



%%
% Smooth the path
switch smoothing 
    case 'on'
        path = getSmoothPath(path, rmin, n_smooth);
end

% Solve for the path
[pathdat,singloc] = qt_bdint_path_solve(path, solseries, []);

% Plot
omega = pathdat.tau - 1i*pathdat.theta;
W = exp(omega);


figure(5); hold all;
plot3(real(pathdat.f),imag(pathdat.f), reim(W)+myeps, '-', 'LineWidth', 2); 
plot3(real(pathdat.f),imag(pathdat.f), 0*reim(W) + minz, '-', 'LineWidth', 2);
xlabel('\phi'); ylabel('\psi');

C = [pathdat.f(:).'; W(:).'];
Cnew = getComplexCurveByArcLength(C, n_arc, reimname);
C = Cnew;
handles = drawPatchForComplexContour(C, patch_width, reimname);
%}

% Routine using meshing --- not as good
%{
base = [-1, 1; 0, 0]/10;
[X,Y,Z] = drawMeshForComplexContour(base, C, 'imag');
h_surf = surf(X,Y,Z);
set(h_surf, 'EdgeColor', 'none');
hold on
for j = 1:size(base, 2)
    h_edge(j) = plot3(X(j,:), Y(j,:), Z(j,:));
    set(h_edge(j), 'Color', 'k', 'LineWidth', 2);
end
%}