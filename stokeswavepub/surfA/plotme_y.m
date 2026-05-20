% Created to make Riemann surface for a spiral
%
% ep = 0.3 using circle.A.A_.ep.0.3
%
% Re-edited Nov. 2015 to get the profile in y which seems more natural

function plotme_y
    
    clear
    F = load('surfA_.ep0.3.mat');
    pathdat = F.pathdat;
    path = F.path;
    
    minz = -0.1;
    dom = [-0.5, 0.5, -0.1, 0.35];
    zeps = 0.005;
    zwidth = 0.01;
%     colorop = 'grey';
    edgecol = 0.7*[1 1 1];
%     basecol = 0.92*[1 1 1];
    basecol = 0.98*[1 1 1];
    outwidth=1;
    
    plottype = 'z';
    reimname = 'imag';
    reim = str2func(reimname);
    
    %%
    close all
    
    figure(1); hold on
    plotplane
    
    % Now run the command
    cd ../../stokeswavecode/bdint
    handles = qt_bdint_path_plotter(path, pathdat, 'PatchWidth', 0.035, 'NumArc', 150, 'ZPlot', reimname, 'ZFunc', plottype);    
    
    % I also want to plot the profile
    f_profile = linspace(-0.5, 0.5, 50);
    [z_profile, ~] = zfree(f_profile, F.solseries);
    
    figure(1); hold on;
    h_profile = plot3(f_profile, 0*f_profile, reim(z_profile), 'k', 'LineWidth', 3);
    
    cd ../../stokeswavepub/surfA
    
    %%
        
%     camproj('perspective')
%     view([-100 20]);

    figure(1)
    
%     delete(handles.surf.edge_right{1})
%     delete(handles.surf.edge_left{1})
    
    %     view([-70 15]);
    view([-20, 20]);
    %     set(gcf, 'Units', 'pixels', 'Position', [60   393   973   555]);
    set(gcf, 'Units', 'pixels', 'Position', [395         140        1395        1051]);
    zlim([minz - zwidth, 0.55])
    xlim([-0.55, 0.55]);
    ylim([-0.1, 0.4]);
    grid on
    xlabel('real potential, \phi', 'FontSize', 24);
    ylabel('imaginary potential, \phi', 'FontSize', 24);
    zlabel('Im(z)', 'FontSize', 24);
    
    
    % Animation for the Manchester talk
    ylim([-0.1, 1e-5]);
    view([0 0]);
    for j = 0:33
        view([-j, j]);
        drawnow
        exportmovie
    end
    
    boo = linspace(1e-5, 0.4, 30);
    for j = 1:length(boo)
        ylim([-0.1, boo(j)]);
        drawnow
        exportmovie
    end       
    
    
    
    figure(99)
    set(gcf, 'Units', 'pixels', 'Position', [60   393   973   250]);
    
    
    function plotplane
        % Bottom plane
        fill3([dom(1) dom(2) dom(2) dom(1)], [dom(3) dom(3) dom(4) dom(4)], (minz - zwidth)*[1 1 1 1], basecol,'LineWidth', outwidth);
        % Left plane
        fill3(dom(1)*[1 1 1 1], [dom(3) dom(3) dom(4) dom(4)], [(minz - zwidth) (minz) (minz) (minz - zwidth)], basecol,'LineWidth', outwidth);
        % Front plane
        fill3([dom(1) dom(2) dom(2) dom(1)], dom(3)*[1 1 1 1], [(minz - zwidth) (minz - zwidth) (minz) (minz)], basecol,'LineWidth', outwidth);
        % Top plane
        fill3([dom(1) dom(2) dom(2) dom(1)], [dom(3) dom(3) dom(4) dom(4)], (minz)*[1 1 1 1], basecol, 'LineWidth', outwidth);
        
        f_A = 1i*0.175;
        plot3([0], imag(f_A)*[1], minz*[1]+zeps, 'ko', 'MarkerFace', 'w', 'MarkerSize', 10, 'LineWidth', 2);
        
        plot3([dom(1), dom(2)], [0, 0], minz*[1 1], 'k--');
        plot3([0, 0], [dom(3), dom(4)], minz*[1 1], 'k--');
        
        plot3(real(pathdat.f),imag(pathdat.f), 0*pathdat.f + minz, 'k-', 'LineWidth', 2);
        
    end
    
end


