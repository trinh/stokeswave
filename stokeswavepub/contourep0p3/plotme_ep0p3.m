function plotme_ep0p3
    % PLOTME code generates contour plots for dualmesh output
    clear
    close all
    
    myshow = 'newfig';
    
    reimname = 'real';
    reim = str2func(reimname);
    
    dirlist = dir('dualmesh*.mat');
    
    exportfig = 'false';
    for j = 1:length(dirlist)
        F = load(dirlist(j).name);
        pathdat = F.pathdat;
        
        plotcontour
        figshift
        title([dirlist(j).name, ', ep = 0.3'], 'Interpreter', 'none');
        
        switch exportfig
            case 'true'
                savefig(gcf, [dirlist(j).name, '.fig'], 'compact');
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
        hold on
        [C, h] = contourf(real(pathdat.down.FF), imag(pathdat.down.FF), reim(omega.down), [-10:0.2:10]);
        %         [C, h] = contourf(real(pathdat.up.FF), imag(pathdat.up.FF), reim(exp(omega.up)), [-10:0.2:10]);
        %         hold on
        %         [C, h] = contourf(real(pathdat.down.FF), imag(pathdat.down.FF), reim(exp(omega.down)), [-10:0.2:10]);
        plot(real(pathdat.test.f), imag(pathdat.test.f), 'LineWidth', 2);
        
        hold off
    end
end