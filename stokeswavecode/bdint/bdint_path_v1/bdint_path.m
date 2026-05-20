clear

% =====================================
%         Load in the solution
F = load('data/bdint_N801_ep0.05.mat');
solbdint = F.solbdint;
F = load('data/series_N100_ep0.05.mat');
solseries = F.solseries;
F = load('data/meshdatbd_ep0.05.mat');
meshdatbd = F.meshdatbd;
% =====================================


% =====================================
%            Path input

% Define the path manually 
% pm = 1i*1e-3;
pm = 1i*1e-3;
path_n1_n1 = [-0.4 + pm, ...
        -0.4 + 1i, ...
        0.4 + 1i, ...
        0.4 + pm, ...
        NaN, ...
        0.4 - pm, ...
        0.4 - 1i, ...
        -0.4 - 1i, ...
        -0.4 - pm];
path_n1_0 = [-0.4 + pm, ...
        -0.4 + 1i, ...
        0.4 + 1i, ...
        0.4 + pm, ...
        NaN, ...
        0.4 - pm, ...
        0.4 - 0.1i, ...
        -0.4 - 0.1i, ...
        -0.4 - pm];
path_n1= [-0.4 + pm, ...
        -0.4 + 1i, ...
        0.4 + 1i, ...
        0.4 + pm];  
path_n0= [-0.4 + pm, ...
        -0.4 + 0.2i, ...
        0.4 + 0.2i, ...
        0.4 + pm];  
path = [path_n1_n1, NaN, path_n1_n1, NaN, path_n1_n1, NaN, path_n1_n1];
% path = path_n1_n1;
path = [path_n1_0, NaN, path_n1_0];
path = [path_n1_n1, NaN, path_n1_n1];
path = [-0.4 + pm, -0.4 + 1i, 0.4 + 1i, 0.4 - 1i, -0.4 - 1i, -0.4 + 1i, ...
    0.4 + 1i, 0.4 - 1i, -0.4 - 1i, -0.4 + 1i];
% path = [path_n1, path_n1, path_n1(1)];
path = path_n0;
     
% Define via mouse
%{
figure(3);
hold on;
view([0, 90]);
path = 0; x_in = []; y_in = []; butt = 1;
while butt == 1
    [x_in(end+1), y_in(end+1), butt] = ginput(1);
end
path = [path, x_in + 1i*y_in];
%}

%path = [0; x_in+1i*y_in];

% =====================================

 
[f, x, y] = bdint_pathsolve(path, solbdint, solseries);
z = x + 1i*y;

%%

figure(1); figshift
% subplot(1, 2, 1); hold all;
plot3(real(f),imag(f),real(x), 'k-', 'LineWidth', 2);
xlabel('\phi'); ylabel('\psi'); zlabel('Re(z)');

% subplot(1, 2, 2); hold all;
% plot3(real(f),imag(f),imag(z), '-', 'LineWidth', 2);
% xlabel('\phi'); ylabel('\psi'); zlabel('Im(z)');

% figure(2); figshift
% subplot(1, 2, 1); hold all;
% plot3(real(f),imag(f),real(y), '-', 'LineWidth', 2);
% xlabel('\phi'); ylabel('\psi'); zlabel('Re(y)');
% 
% subplot(1, 2, 2); hold all;
% plot3(real(f),imag(f),imag(y), '-', 'LineWidth', 2);
% xlabel('\phi'); ylabel('\psi'); zlabel('Im(y)');

% figure(3); figshift
% subplot(1, 2, 1); hold all;
% plot3(real(f),imag(f),real(z), '-', 'LineWidth', 2);
% xlabel('\phi'); ylabel('\psi'); zlabel('Re(z)');
% 
% subplot(1, 2, 2); hold all;
% plot3(real(f),imag(f),imag(z), '-', 'LineWidth', 2);
% xlabel('\phi'); ylabel('\psi'); zlabel('Im(z)');
% 
