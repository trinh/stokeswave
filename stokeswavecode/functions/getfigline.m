function [x, y, lineH] = getfigline
% GETFIGLINE gets the line data from the current figure and returns the x
% data, y data and the line figure handle
%
% Then access via x{j}, y{j}, lineH{j}

h = gcf;

lineH = findobj(h, 'type', 'line');
x = get(lineH, 'Xdata');
y = get(lineH, 'Ydata');
