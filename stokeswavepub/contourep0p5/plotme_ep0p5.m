function plotme_ep0p5
    % PLOTME code generates contour plots for dualmesh output
    clear
    close all
    
    myshow = 'newfig';
    
    reimname = 'real';
    reim = str2func(reimname);
    
    dirlist = dir('./100x200/dualmesh*.mat');
    
    exportfig = 'false';
    exportfile = 'true';
    for j = 1:length(dirlist)
        F = load(['./100x200/', dirlist(j).name]);
        [pathstr,name,ext] = fileparts(dirlist(j).name);
        
        pathdat = F.pathdat;
        
        
        
%                 plotcontour
        plotpub
        figshift
        
        
        switch exportfig
            case 'true'
                savefig(gcf, [name, '.fig'], 'compact');
        end
        switch exportfile
            case 'true'
                export_fig([name, '.pdf']);
        end
        
    end
    
    function plotcontour
        switch myshow
            case 'newfig'
                figure
            case 'subfig'
        end
        omega.up = pathdat.up.taumat - 1i*pathdat.up.thetamat;
        omega.down = pathdat.down.taumat - 1i*pathdat.down.thetamat;
        [C, h] = contourf(real(pathdat.up.FF), imag(pathdat.up.FF), reim(omega.up), [-10:0.2:10]);
        %         set(h, 'EdgeColor', 0.5*[1 1 1]);
        hold on
        [C, h] = contourf(real(pathdat.down.FF), imag(pathdat.down.FF), reim(omega.down), [-10:0.2:10]);
        %         set(h, 'EdgeColor', 0.5*[1 1 1]);
        %         [C, h] = contourf(real(pathdat.up.FF), imag(pathdat.up.FF), reim(exp(omega.up)), [-10:0.2:10]);
        %         hold on
        %         [C, h] = contourf(real(pathdat.down.FF), imag(pathdat.down.FF), reim(exp(omega.down)), [-10:0.2:10]);
        plot(real(pathdat.test.f), imag(pathdat.test.f), 'LineWidth', 2);
        hold off
        title(name, 'Interpreter', 'none');
    end
    function plotpub
        figure
        omega.up = pathdat.up.taumat - 1i*pathdat.up.thetamat;
        omega.down = pathdat.down.taumat - 1i*pathdat.down.thetamat;
        
       
        [C, h] = contour(real(pathdat.up.FF), imag(pathdat.up.FF), reim(omega.up), [-10:0.2:10]);
        set(h, 'EdgeColor', 0.3*[1 1 1], 'LineWidth', 2);
        hold on
        [C, h] = contour(real(pathdat.down.FF), imag(pathdat.down.FF), reim(omega.down), [-10:0.2:10]);
        set(h, 'EdgeColor', 0.3*[1 1 1], 'LineWidth', 2);
        hold off
        ylim([-1, 0.5]);
        set(gca, 'Visible', 'off');
        
    end
end