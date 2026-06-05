%A script to stop export_fig bug with darkmode

set(gcf,'Color','w')
set(gca,'Color','w')

set(gca,'XColor','k','YColor','k','ZColor','k') % axes lines + ticks
set(findall(gcf,'Type','text'),'Color','k')     % all text
set(findall(gcf,'Type','axes'),'Color','w')     % all axes backgrounds

set(gcf,'InvertHardcopy','off')

lines = findall(gcf,'Type','line');

for k = 1:length(lines)
    c = get(lines(k),'Color');
    
    % MATLAB default blue ≈ [0 0.4470 0.7410]
    if norm(c - [0 0.4470 0.7410]) < 1e-3
        set(lines(k),'Color','k')
    end
end