% Test to extract sam data

clear
close all

open ../singularitylocations_withseries.fig
[x, y, lineH] = getfigline;

% First on is the blue curve
ep = x{1};
loc = y{1};

figure
plot(log(1 - ep), log(loc), 'bo');
hold on
plot(log(1 - ep), 1.1082*log(1 - ep) - 1.8, 'k--');
xlabel('log(1 - ep)');
ylabel('location');
title('Trend as epsilon to one');

%%

clear
close all

open ./first_sing_loc.fig
[x, y, lineH] = getfigline;

% First on is the blue curve
ep = x;
loc = y;

figure
plot(log(ep), loc, 'bo');
xlabel('log(epsilon)')
ylabel('location');
hold on
plot(log(ep), -0.15672*log(ep) -0.013747, 'r');
title('Trend as epsilon to zero');

