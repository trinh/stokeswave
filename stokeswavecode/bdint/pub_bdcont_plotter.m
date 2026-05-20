% Code to generate analytic continuation paths along arbitrary contours
function pub_bdcont_plotter
clear
close all

minz = -8;
dom = [-1, 1, -1, 1];
zeps = 0.05;
zwidth = 0.5;
colorop = 'grey';
edgecol = 0.7*[1 1 1];
basecol = 0.92*[1 1 1];
outwidth=2;

% =====================================
%         Load in the solution
F = load('../../stokeswavedata/series_N50_ep0.1.mat');
solseries = F.solseries;
% =====================================


%% Universal constants
smoothing = 'on';
reimname = 'real';
reim = str2func(reimname);
figsize = [600   802   842   543];
figsize = [600   570  1000 750];



%% Parameters
% =====================================
%            Path input
path = get_path('pub_AA_');
rmin = 0.1;
n_smooth = 3;
n_arc = 100;
patch_width = 0.12;


%%
% Smooth the path
% switch smoothing
%     case 'on'
%         path = getSmoothPath(path, rmin, n_smooth);
% end

% Solve for the path
[pathdat,singloc] = qt_bdint_path_solve(path, solseries, []);

% Plot
omega = pathdat.tau - 1i*pathdat.theta;
W = exp(omega);

% Load in data into a 2xm matrix
C = [pathdat.f(:).'; W(:).'];


figure(5); figshift; hold on

plotplane % Plot the plane
% plotdash  % Plot the dashes
plotsurf  % Plot the surface

hold on
plot3(real(pathdat.f),imag(pathdat.f), reim(W)+zeps, 'k-', 'LineWidth', 1); 

% Make graph pretty
set(gcf, 'Units', 'pixels', 'Position', figsize);
set(gca,'ztick',[-8:4:8]);
set(gca, 'FontSize', 18);

% xlabel('\phi'); ylabel('\psi');
zlim([-8-zwidth, 8]);
view([-50, 40]);
grid on
shg

% One plot to plot the curve as a function of arclength
% [Cnew,s] = getComplexCurveByArcLength(C, n_arc, reimname);
% figure(6)
% plot(s, reim(Cnew(2,:)));



    function plotsurf
        Cnew = getComplexCurveByArcLength(C, n_arc, reimname, 'flatten');
        handles = drawPatchForComplexContour(Cnew, patch_width, reimname);
        for j = 1:length(handles.edge_right)
            delete(handles.edge_right{j});
            delete(handles.edge_left{j});
        end
    end


    function plotdash
        dashnum = 100;
        Cflat = getComplexCurveByArcLength(C, dashnum, reimname, 'flatten');
%         for j = 1:size(Cflat, 2)
% %             plot3(real(Cflat(1,j))*[1 1], imag(Cflat(1,j))*[1 1], [minz, reim(Cflat(2,j))], 'k--', 'Color', 0.5*[1 1 1], 'LineWidth', 0.5);
%                 plot3(real(Cflat(1,j))*[1 1], imag(Cflat(1,j))*[1 1], [minz, reim(Cflat(2,j))], 'k', 'Color', 0.7*[1 1 1], 'LineWidth', 0.5);
%         end
        for j = 1:size(Cflat,2) - 1
            x = real([Cflat(1,j) Cflat(1,j) Cflat(1,j+1) Cflat(1,j+1)]);
            y = imag([Cflat(1,j) Cflat(1,j) Cflat(1,j+1) Cflat(1,j+1)]);
            z = [minz+zeps reim(Cflat(2,j)) reim(Cflat(2,j+1)) minz+zeps];
            p = patch(x,y,z,0.5*[1 1 1]);
            set(p, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
%             set(p, 'FaceAlpha', 0.5);
        end
    end

    function plotplane
        % Bottom plane
        fill3([dom(1) dom(2) dom(2) dom(1)], [dom(3) dom(3) dom(4) dom(4)], (minz - zwidth)*[1 1 1 1], basecol,'LineWidth', outwidth);
        % Top plane
        fill3([dom(1) dom(2) dom(2) dom(1)], [dom(3) dom(3) dom(4) dom(4)], (minz)*[1 1 1 1], basecol, 'LineWidth', outwidth);
        % Left plane
        fill3(dom(1)*[1 1 1 1], [dom(3) dom(3) dom(4) dom(4)], [(minz - zwidth) (minz) (minz) (minz - zwidth)], basecol,'LineWidth', outwidth);
        % Front plane
        fill3([dom(1) dom(2) dom(2) dom(1)], dom(3)*[1 1 1 1], [(minz - zwidth) (minz - zwidth) (minz) (minz)], basecol,'LineWidth', outwidth);
        %         plot3([-2 -1 0], [0 0 0], minz*[1 1 1], 'ko', 'MarkerFace', 'w', 'MarkerSize', 10, 'LineWidth', 2);
        
        plot3([dom(1), dom(2)], [0, 0], minz*[1 1], 'k--');
        
        plot3(real(pathdat.f),imag(pathdat.f), 0*reim(W) + minz, 'k-', 'LineWidth', 2);
        
    end
end
