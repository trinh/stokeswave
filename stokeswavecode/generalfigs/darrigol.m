% Several wave profiles plotted s.t. their singularities overlap.

clear
close all

N = 80; % Truncation
M = 8; % Epsilon iterations (7 for figure)
K = 200; % Plotting phi mesh
kap = 0;
grav_coeff = 1;

eps_sq = linspace(0.1,0.8, M);
%eps_sq = linspace(0.6, 0.8, M);
%eps_sq = [0.1, 0.15, 0.2, 0.3, 0.4];
eps = eps_sq.^(0.5);

p=0.1;
q = 0.4;
path = [p, p+q*1i, -p+q*1i, -p, p, p+q*1i, -p+q*1i, -p, p];

% Matrix for z data
zmat = [];

figure(1); hold all
xlabel('x');
ylabel('y');
for k=1:length(eps_sq)
% Initial guess - last two are mu and B (in that order!)
x0 = zeros(1,N-1);
x0(N) = 1;
x0(N+1) = 0.1;
ep = eps_sq(k)^(0.5);

% Solve the ODE
options = optimset('Display', 'iter');
fwd = @(x)series_bern(x, ep, N, kap, grav_coeff);
[est,fval] = fsolve(fwd,x0,options);

% Save all the data
solseries = [];
solseries.N = N;
solseries.an = est(1:N-1);
solseries.mu = est(N);
solseries.B = est(N+1);
solseries.ep = ep;
solseries.kap = kap;
solseries.ycrest = sum(est(1:N-1)./(2*pi*(1:N-1)));

% Calculate the wave profile
phi = linspace(-0.5, 0.5, K);
[z,~] = zfree(phi, solseries);

[f,g,gp,~,singloc] = gbern_contour(path,solseries);
path_sing = [p, singloc];
[f2,g2,gp2,~,~] = gbern_contour(path_sing,solseries);

% Then we calculate the y value at the singularity. g(f) = z(f) + z(-f)
[z_sing_minus,~] = zfree(-singloc, solseries);
z_sing(k) = g2(end) - z_sing_minus;

location(k) = singloc;
B(k) = solseries.B;
mu(k) = solseries.mu;
H(k) = imag(-zfree(0.5, solseries)+zfree(0, solseries));

% Phil modification to plot two periods
x_add = [real(z), real(z) + 1];
y_add = [imag(z), imag(z)];

%plot3(phi,zeros(1, length(phi)),imag(z)); % 3D version
% plot(real(z),imag(z)); % Unshifted
plot(x_add, y_add - imag(z_sing(k)));

zmat{end+1} = z;
end


%plot(real(z_sing),imag(z_sing),'b*','MarkerSize', 20); % Unshifted version
%plot(real(location),imag(location),'b*','MarkerSize', 20); % 3D Version

% data.B = B;
% data.mu = mu;
% data.H = H;
% data.eps_sq = eps_sq;
% data.location_psi = location;
% data.location_y = imag(z_sing);
% 
% save('darrigoldat.mat', 'data');

% save darrigoldat.mat;