% Plot the upper mesh for the publication using regular series
function meshupper_plotpub

close all

dom = [-0.5, 0.50, -1, 1];
minz = -1;
zeps = 0.05;
zwidth=0.03;
colorop = 'grey';
edgecol = 0.7*[1 1 1];
basecol = 0.92*[1 1 1];

figure(1); clf(1); hold on
outwidth = 2;

% load in the lower plane data
F = load('../../../stokeswavedata/recdat_rec4_ep0.1.mat');
recdat = F.recdat;

f_d = recdat(1).f;
z_d = recdat(1).z;
f_u = recdat(2).f;
z_u = recdat(2).z;

% Trim
z_d(imag(f_d) > -1e-2) = NaN;
f_d(imag(f_d) > -1e-2) = NaN;
z_u(imag(f_u) < 1e-2) = NaN;
f_u(imag(f_u) < 1e-2) = NaN;

% load in surface data
F = load('../../../stokeswavedata/meshdat_M100_ep0.1.mat');
meshdat = F.meshdat;
% trim the surface data around the cut
UHP = meshdat.UHP(1);
f_left = UHP.f; z_left = UHP.z;
f_right = UHP.f; z_right = UHP.z;
z_left(real(f_left) > 0) = NaN;
z_right(real(f_right) < 0) = NaN;
LHP = meshdat.LHP;

% ------------- PLOTTING COMMANDS -------------------

% plot the lower plane
% plotplane

% plot the surface
% plotsurf

% plot the surface contours
plotcont

grid on
view([-38, 24]);
set(gcf, 'Position', [500 500 600 700]);
set(gca,'xtick',[-0.5:0.25:0.5])
set(gca, 'Visible', 'off');
xlim([-0.5, 0.5]);
ylim([-1, 1]);
zlim([-1.2, 0.6]);

% ------------- END PLOT -------------------

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
        
        
        
        plot3(real(f_d), imag(f_d), 0*real(z_d)+minz, 'k', 'LineWidth', 4);
        plot3(real(f_u), imag(f_u), 0*real(z_u)+minz, 'k', 'LineWidth', 4);
    end

    function plotsurf
        [h_left, hwire_left] = plotme(f_left, real(z_left), 'interp', 10); hold all
        [h_right, hwire_right] = plotme(f_right, real(z_right), 'interp', 10);
        [h_lower, hwire_lower] = plotme(LHP.f, real(LHP.z), 'interp', 10);
        
        edge_color = 0.7*[1 1 1];
        
        set(h_left, 'FaceColor', 'w');
        set(h_right, 'FaceColor', 'w');
        set(h_lower, 'FaceColor', 'w');
        set(hwire_left, 'EdgeColor', edge_color);
        set(hwire_right, 'EdgeColor', edge_color);
        set(hwire_lower, 'EdgeColor', edge_color);
    end

    function plotcont
        plot3(real(f_d), imag(f_d), real(z_d), 'k', 'LineWidth', 4);
        plot3(real(f_u), imag(f_u), real(z_u), 'k', 'LineWidth', 4);
    end
end