clear
close all

K = 250; % Phi mesh
L = 50; % z(f) integrate mesh

% I think it's easier to load in data for pre-made figures rather than making on
% the go. Singularity location data from singularitylocation_highersheets
% with steep rectangle data. 

% To do the shifted (rather than 3D) need to take the z values from
% zpath_rectangle_steep.m (by inspecting recdat and interpolating). This is not very precise
% but it is possible to make it precise by adapting zpath code.

% Matrix to hold z data
zmat = [];

figure(1); hold all

% Do the Non-Steep waves first (singularity data taken from darrigol.m)

% ========== ep=0.6 ========== %
F = load('../../stokeswavedata/series_N100_ep0.6.mat');
solseries = F.solseries;
phi = linspace(-0.5,0.5,K);
z = zfree(phi, solseries);
singloc(1) = 0.0672;

plot3(phi, zeros(1,length(phi)), imag(z)); % 3D Version

zmat{end+1} = z;

% ========== ep=0.7 ========== %
F = load('../../stokeswavedata/series_N100_ep0.7.mat');
solseries = F.solseries;
phi = linspace(-0.5,0.5,K);
z = zfree(phi, solseries);
singloc(2) = 0.0440;

plot3(phi, zeros(1,length(phi)), imag(z)); % 3D Version

zmat{end+1} = z;

% Now the steep ones.

% ========== ep=0.8 ========== %
F = load('../../stokeswavedata/series_steep_N100_ep0.8.mat');
steepsolseries = F.steepsolseries;
phi = linspace(-0.5,0.5,K);
z = zfree_steep(phi, steepsolseries, L);
singloc(3) = 0.0248;

plot3(phi, zeros(1,length(phi)), imag(z)); % 3D Version

zmat{end+1} = z;

% ========== ep=0.9 ========== %
F = load('../../stokeswavedata/series_steep_N100_ep0.9.mat');
steepsolseries = F.steepsolseries;
phi = linspace(-0.5,0.5,K);
z = zfree_steep(phi, steepsolseries, L);
singloc(4) = 0.0096;

plot3(phi, zeros(1,length(phi)), imag(z)); % 3D Version

zmat{end+1} = z;

% ========== ep=0.95 ========== %
F = load('../../stokeswavedata/series_steep_N100_ep0.95.mat');
steepsolseries = F.steepsolseries;
phi = linspace(-0.5,0.5,K);
z = zfree_steep(phi, steepsolseries, L);
singloc(5) = 0.0040;

plot3(phi, zeros(1,length(phi)), imag(z)); % 3D Version

zmat{end+1} = z;

% ================================ %

plot(zeros(1,length(singloc)),real(singloc),'b*','MarkerSize', 20); % 3D Version

save darrigoldat_steep.mat;