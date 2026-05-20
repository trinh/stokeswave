function h = plot_arcprofile(f, u, path)
% PLOT_ARCPROFILE plots u(f) as a function of the arclength measured via
% the complex curve f. In addition, given a path, it plots dashed lines
% indicating where each node of the path lies

% labels = 0; 

figure(99); clf(99);
set(gcf, 'Units', 'pixels', 'Position',[832 435 1000 300]);
arcl = get_2dcumarc(f);

%{
% Remove NaN from path
path_f = path.f;
path_f(isnan(path_f) == 1) = [];
arcl_path = get_2dcumarc(path_f);

subplot(1, 4, 1)
% We want to interpolate for points along the curve
% But first remove elements that have zero arclength separation
tmp_arcl = arcl; tmp_arcl(diff(arcl) == 0) = [];
tmp_u = u; tmp_u(diff(arcl) == 0) = [];


uinterp = interp1(tmp_arcl, tmp_u, arcl_path);
plot3(real(path_f), imag(path_f), uinterp);
%}


% if labels == 1
%     for j = 1:length(path_f)
%         if mod(j, 5) == 0
%             text(real(path_f(j)), imag(path_f(j)), uinterp(j), num2str(j), 'FontSize', 16, 'VerticalAlignment', 'bottom');
%         elseif mod(j, 5) == 1
%             text(real(path_f(j)), imag(path_f(j)), uinterp(j), num2str(j), 'FontSize', 16, 'VerticalAlignment', 'top');
%         elseif mod(j, 5) == 2
%             text(real(path_f(j)), imag(path_f(j)), uinterp(j), num2str(j), 'FontSize', 16, 'HorizontalAlignment', 'right');
%         elseif mod(j, 5) == 3
%             text(real(path_f(j)), imag(path_f(j)), uinterp(j), num2str(j), 'FontSize', 16, 'HorizontalAlignment', 'left');        
%         elseif mod(j, 5) == 4
%             text(real(path_f(j)), imag(path_f(j)), uinterp(j), num2str(j), 'FontSize', 16, 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right');
%         end    
%     end
% end

% subplot(1, 4, [2:4]);
h = [];
h(1) = plot(arcl, real(u), 'k'); 
hold on
h(2) = plot(arcl, imag(u), 'k--');

% Mark the location of singstart
% h(3) = plot(arcl(path.singstart)*[1, 1], [min([real(u), imag(u)]), max([real(u), imag(u)])], 'k--');

% if labels == 1
%     hold on
%     for j = 1:length(arcl_path)
%         h(end+1) = plot(arcl_path(j)*[1, 1], [min(u), max(u)], 'k--');    
%         text(arcl_path(j), max(u), num2str(j), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'left');    
%     end
%     hold off
% end

title('u(arclength)');
% keyboard