% Plot Stokes wave for pub
% Run darrigol_steep.m
clear
close all

% load 'darrigoldat_steep.mat';
load 'darrigoldat.mat';

figure(1); clf(1); hold on;
for j = 1:length(zmat)
    z = zmat{j};

    x_add = [real(z), real(z) + 1];
    y_add = [imag(z), imag(z)];

    plot(x_add, y_add - imag(z_sing(j)), 'k');
end

%%
clear
load 'darrigoldat_steep_sq.mat'
for j = 1:length(zmat)
    z = zmat{j};

    x_add = [real(z), real(z) + 1];
    y_add = [imag(z), imag(z)];

    plot(x_add, y_add - ysingloc(j), 'k--');
end

%% 

clear
close all

% load 'darrigoldat_steep.mat';
load 'darrigoldat.mat';

figure(1); clf(1); hold on;
for j = 1:length(zmat)
    z = zmat{j};

    x_add = [real(z), real(z) + 1];
    y_add = [imag(z), imag(z)];

    plot(x_add, y_add - imag(z(1)), 'k');
    plot(0, imag(z_sing(j) - z(1)), 'o');
end

%%

clear
load 'darrigoldat_steep_sq.mat'
for j = 1:1
    z = zmat{j};

    x_add = [real(z), real(z) + 1];
    y_add = [imag(z), imag(z)];

    plot(x_add, y_add - imag(z(1)), 'k--'); 
    plot(0, ysingloc(j) - imag(z(1)), 'o');
end