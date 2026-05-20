clear
close all

K = 500; % Phi mesh
L = 1000; % z(f) integrate mesh
phi = linspace(-0.5,0.5,K);

% Matrix to hold z data
zmat = [];
% Note these are specified
eps_sq_mat = [0.9, 0.94, 0.99];

figure(1); hold all

% Steep wave data in terms of ep^2.
% Singularity location data comes from
% singularitylocation/singularitylocation_highersheets
% ysingloc is found from inspecting rectangular data in
% zpath_rectangle_steep.m

% ========== ep=0.9487 (ep^2 = 0.9)========== %
F = load('../../stokeswavedata/series_steep_N100_ep0.9487.mat');
steepsolseries = F.steepsolseries;
mu(1) = steepsolseries.mu;
H(1) = imag(-zfree_steep(0.5, steepsolseries, L)+zfree_steep(0, steepsolseries, L));
z = zfree_steep(phi, steepsolseries, L);
% singloc(1) = 0.0035; % Sam old figure
% singloc(1) = 0.000060289882+0.003493677829i; % Phil correction --- Re(f) off by 1e-5
singloc(1) = -0.000000414956+0.000005065235i; % Phil correction v2 using more accurate recdat---wow!
% Even if the singularity is off, use only the imaginary part
ysingloc(1) = imag(zfree_steep(imag(singloc(1))*1i, steepsolseries, L));

%plot3(phi, zeros(1,length(phi)), imag(z)); % 3D Version
plot(real(z), imag(z-ysingloc(1)));

zmat{end+1} = z;

% ========== ep=0.9695 (ep^2 = 0.94) ========== %
F = load('../../stokeswavedata/series_steep_N100_ep0.9695.mat');
steepsolseries = F.steepsolseries;
mu(2) = steepsolseries.mu;
H(2) = imag(-zfree_steep(0.5, steepsolseries, L)+zfree_steep(0, steepsolseries, L));
z = zfree_steep(phi, steepsolseries, L);
% singloc(2) = 0.0021;
% singloc(2) = 0.000343100585+0.002117919166i; % Phil correction, off by 1e-4
singloc(2) = 0.000326345996+0.002039919745i; % Phil correction v2
% ysingloc(2) = 0.013;
ysingloc(2) = imag(zfree_steep(imag(singloc(2))*1i, steepsolseries, L));
%plot3(phi, zeros(1,length(phi)), imag(z)); % 3D Version
plot(real(z), imag(z-ysingloc(2)));

zmat{end+1} = z;

% ========== ep=0.995 (ep^2 = 0.99) ========== %
F = load('../../stokeswavedata/series_steep_N100_ep0.995.mat');
steepsolseries = F.steepsolseries;
mu(3) = steepsolseries.mu;
H(3) = imag(-zfree_steep(0.5, steepsolseries, L)+zfree_steep(0, steepsolseries, L));
z = zfree_steep(phi, steepsolseries, L);
singloc(3) = 0.013312050205+0.000404408608i; % Way off!!!
singloc(3) = 0.000231581980+0.000112374823i; % Phil correction v2
ysingloc(3) = imag(zfree_steep(imag(singloc(3))*1i, steepsolseries, L));

%plot3(phi, zeros(1,length(phi)), imag(z)); % 3D Version

zmat{end+1} = z;

% ================================ %

%plot(zeros(1,length(singloc)),real(singloc),'b*','MarkerSize', 20); % 3D Version

data.mu = mu;
data.H = H;
data.eps_sq = [0.9, 0.94, 0.99];
data.location_psi = singloc;
data.location_y = ysingloc;

save 'darrigoldat_steep_sq.mat'
