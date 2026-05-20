% Created to make Riemann surface for a spiral
%
% ep = 0.3 using circle.A.A_.ep.0.3

function plotme
    
    clear
    F = load('surf.A.A_.ep0.3_smooth2.mat');
    pathdat = F.pathdat;
    path = F.path;
    
    minz = -0.3;
    dom = [-0.2, 0.2, -0.35, 0.35];
    zeps = 0.02;
    zwidth = 0.08;
    colorop = 'grey';
    edgecol = 0.7*[1 1 1];
    basecol = 0.92*[1 1 1];
    outwidth=1;
    
    reimname = 'real';
    reim = str2func(reimname);
    
    %%
    close all
    
    figure(1); hold on
    plotplane
    
    % Now run the command
    cd ../../stokeswavecode/bdint
    handles = qt_bdint_path_plotter(path, pathdat, 'PatchWidth', 0.02, 'NumArc', 150, 'ZPlot', 'imag', 'ZFunc', 'omega');
    cd ../../stokeswavepub/surfA_
    
    %%
        
%     camproj('perspective')
%     view([-100 20]);

    figure(1)
    
%     delete(handles.surf.edge_right{1})
%     delete(handles.surf.edge_left{1})
    
    view([100 50]);
    set(gcf, 'Units', 'pixels', 'Position', [60   393   973   555]);
    zlim([minz - zwidth, 1.6])
    xlim([-0.2, 0.2]);
    ylim([-0.35, 0.35]);
    
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
        plot3([0, 0], imag(f_A)*[1 -1], minz*[1 1]+zeps, 'ko', 'MarkerFace', 'w', 'MarkerSize', 10, 'LineWidth', 2);
        
        plot3([dom(1), dom(2)], [0, 0], minz*[1 1], 'k--');
        plot3([0, 0], [dom(3), dom(4)], minz*[1 1], 'k--');
        
        plot3(real(pathdat.f),imag(pathdat.f), 0*pathdat.f + minz, 'k-', 'LineWidth', 2);
        
    end
    
end


