function plotme
    F = load('meshAA_ep0.5_50100.mat');
    G = load('pathAA_.B_ep0.5.mat');
    meshpathdat = F.pathdat;
    pathdat = G.pathdat;
    path = G.path;
    
    
    minz = -3.5;
    dom = [-1.2, 1.2, -0.8, 0.6];
    zeps = 0.08;
    zwidth = 0.1;
    colorop = 'grey';
    edgecol = 0.7*[1 1 1];
    meshedgecol = 0.35*[1 1 1];
    basecol = 0.92*[1 1 1];
    outwidth=1;
    
    zmeshmin = -2;
    zmeshmax = 2.8;
    
    reimname = 'imag';
    reim = str2func(reimname);
    
    
    
    close all
    figure(1); clf(1); hold on
    set(gcf, 'Units', 'pixels', 'Position', [620 268 1457 1070]);
    zlim([minz - zwidth, 3])
    xlim([-1.2, 1.2]);
    ylim([-0.8, 0.6]);
    view([-10, 20]);
    grid on
%     colormap gray
    
    plotplane
%     grid off
%     set(gca, 'Visible', 'off');

    plotsurf_up;
%     w = wireframe(X,Y,Z,4, .5, meshedgecol);

    plotsurf_down
%     w = wireframe(X,Y,Z,4, .5, meshedgecol);
%     brighten(0.5);
    
    
    % Now run the command
    cd ../../stokeswavecode/bdint
    handles = qt_bdint_path_plotter(path, pathdat, 'PatchWidth', 0.07, 'NumArc', 180, 'ZPlot', 'imag', 'ZFunc', 'dfdz');
    
    delete(handles.surf.edge_right{1})
    delete(handles.surf.edge_left{1})
    cd ../../stokeswavepub/surfAA_.B_/
    
    function plotsurf_up
        X = real(meshpathdat.up.FF);
        Y = imag(meshpathdat.up.FF);
        Z = reim(exp(meshpathdat.up.taumat - 1i*meshpathdat.up.thetamat)) - zeps;
        
        ind = find(Y(:,1) > dom(4));
        
        X(ind,:) = [];
        Z(ind,:) = [];
        Y(ind,:) = [];
        
        p = surf(X, Y, Z);
%         set(p, 'EdgeColor', 'none', 'FaceColor', 0.9*[1 1 1]);        
        set(p, 'EdgeColor', 0.8*[1 1 1], 'FaceColor', 0.9*[1 1 1]);        
        set(p, 'EdgeColor', 0.5*[1 1 1], 'FaceColor', 0.9*[1 1 1]);        
    end
    function plotsurf_down
        X = real(meshpathdat.down.FF);
        Y = imag(meshpathdat.down.FF);
        Z = reim(exp(meshpathdat.down.taumat - 1i*meshpathdat.down.thetamat)) - zeps;

        ind = find(Y(:,1) < dom(3));

        X(ind,:) = [];
        Z(ind,:) = [];
        Y(ind,:) = [];

        Z(find(Z < zmeshmin)) = zmeshmin;
        Z(find(Z > zmeshmax)) = zmeshmax;

        p = surf(X, Y, Z);
%         set(p, 'EdgeColor', 'none', 'FaceColor', 0.9*[1 1 1]);
set(p, 'EdgeColor', 0.8*[1 1 1], 'FaceColor', 0.9*[1 1 1]);  
set(p, 'EdgeColor', 0.5*[1 1 1], 'FaceColor', 0.9*[1 1 1]);        
    end
    function plotplane
        % Bottom plane
        fill3([dom(1) dom(2) dom(2) dom(1)], [dom(3) dom(3) dom(4) dom(4)], (minz - zwidth)*[1 1 1 1], basecol,'LineWidth', outwidth);
        % Left plane
        fill3(dom(1)*[1 1 1 1], [dom(3) dom(3) dom(4) dom(4)], [(minz - zwidth) (minz) (minz) (minz - zwidth)], basecol,'LineWidth', outwidth);
        % Right plane
        fill3(dom(2)*[1 1 1 1], [dom(3) dom(3) dom(4) dom(4)], [(minz - zwidth) (minz) (minz) (minz - zwidth)], basecol,'LineWidth', outwidth);
        % Front plane
        fill3([dom(1) dom(2) dom(2) dom(1)], dom(3)*[1 1 1 1], [(minz - zwidth) (minz - zwidth) (minz) (minz)], basecol,'LineWidth', outwidth);
        % Top plane
        fill3([dom(1) dom(2) dom(2) dom(1)], [dom(3) dom(3) dom(4) dom(4)], (minz)*[1 1 1 1], basecol, 'LineWidth', outwidth);
        
        f_A = 1i*0.175;
        plot3([0, 0, -1 1], imag(f_A)*[1 -1, 1, 1], minz*[1 1, 1, 1]+zeps, 'ko', 'MarkerFace', 'w', 'MarkerSize', 10, 'LineWidth', 2);
        
        f_B = 0.545 - 0.387i;
        plot3(real(f_B)*[-1 1], imag(f_B)*[1 1], minz*[1 1]+zeps, 'ko', 'MarkerFace', 'w', 'MarkerSize', 10, 'LineWidth', 2);
        
        plot3([dom(1), dom(2)], [0, 0], minz*[1 1], 'k--');
        plot3([0, 0], [dom(3), dom(4)], minz*[1 1], 'k--');
        
        plot3([1, 1], [dom(3), dom(4)], minz*[1 1], 'k--');
        plot3([-1, -1], [dom(3), dom(4)], minz*[1 1], 'k--');
        
        plot3([0.5, 0.5], [dom(3), dom(4)], minz*[1 1], 'k-.');
        plot3([-0.5, -0.5], [dom(3), dom(4)], minz*[1 1], 'k-.');
        
        plot3(real(pathdat.f),imag(pathdat.f), 0*pathdat.f + minz, 'k-', 'LineWidth', 2);
        
    end
    
end