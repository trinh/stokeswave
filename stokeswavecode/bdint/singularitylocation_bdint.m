% Finding the distance from the origin of various singularities as a
% function of epsilon

%clear
close all

N = 50; % Truncation
M = 4; % Epsilon iterations
kap = 0;
grav_coeff = 1;

eps = linspace(0.3, 0.7, M);

% I should write this in a loop - getting a bit ridiculous now


% ========= FIRST SINGULARITY ================%
path = get_path('first_singularity');
sing_start = 1;

for k=1:M
% Initial guess - last two are mu and B (in that order!)
x0 = zeros(1,N-1);
x0(N) = 1;
x0(N+1) = 0.1;
ep = eps(k);
est = zeros(1,N+1);

% Solve the ODE
%mytol = 1e-13;
%options = optimset('Display', 'iter', 'TolFun', mytol, 'TolX', mytol);
options = optimset('Display', 'iter');
fwd = @(x)series_bern(x, ep, N, kap, grav_coeff);
[est,fval] = fsolve(fwd,x0,options); % Call solver

% Save all the data
solseries = [];
solseries.N = N;
solseries.an = est(1:N-1);
solseries.mu = est(N);
solseries.B = est(N+1);
solseries.ep = ep;
solseries.kap = kap;
solseries.ycrest = sum(est(1:N-1)./(2*pi*(1:N-1)));

[~, singloc(1,k)] = qt_bdint_path_solve_sam_singloc(path, solseries, sing_start); % Gets the singularity location
%[~, singcoeff(1,k)] = qt_bdint_path_solve_sam_singcoeff(path, solseries, sing_start, -singloc(k)); % Gets the path and the singularity coefficient
singloc(1,k) = -singloc(1,k);

end


%======== SECOND SINGULARITY ========$

path = get_path('nonphys_singularity');
sing_start = 4;

for k=1:M
% Initial guess - last two are mu and B (in that order!)
x0 = zeros(1,N-1);
x0(N) = 1;
x0(N+1) = 0.1;
ep = eps(k);
est = zeros(1,N+1);

% Solve the ODE
mytol = 1e-13;
options = optimset('Display', 'iter', 'TolFun', mytol, 'TolX', mytol);
fwd = @(x)series_bern(x, ep, N, kap, grav_coeff);
[est,fval] = fsolve(fwd,x0,options); % Call solver

% Save all the data
solseries = [];
solseries.N = N;
solseries.an = est(1:N-1);
solseries.mu = est(N);
solseries.B = est(N+1);
solseries.ep = ep;
solseries.kap = kap;
solseries.ycrest = sum(est(1:N-1)./(2*pi*(1:N-1)));

[~, singloc(2,k)] = qt_bdint_path_solve_sam_singloc(path, solseries, sing_start); % Gets the singularity location
%[~, singcoeff(2,k)] = qt_bdint_path_solve_sam_singcoeff(path, solseries, sing_start, -singloc(k)); % Gets the path and the singularity coefficient

singloc(2,k) = -singloc(2,k);
end



%======== LEFT DIAGONAL ========%

path = get_path('left_diag_singularity');
sing_start = 8;

for k=1:M
% Initial guess - last two are mu and B (in that order!)
x0 = zeros(1,N-1);
x0(N) = 1;
x0(N+1) = 0.1;
ep = eps(k);
est = zeros(1,N+1);

% Solve the ODE
%mytol = 1e-13;
%options = optimset('Display', 'iter', 'TolFun', mytol, 'TolX', mytol);
options = optimset('Display', 'iter');
fwd = @(x)series_bern(x, ep, N, kap, grav_coeff);
[est,fval] = fsolve(fwd,x0,options); % Call solver

% Save all the data
solseries = [];
solseries.N = N;
solseries.an = est(1:N-1);
solseries.mu = est(N);
solseries.B = est(N+1);
solseries.ep = ep;
solseries.kap = kap;
solseries.ycrest = sum(est(1:N-1)./(2*pi*(1:N-1)));

[~, singloc(3,k)] = qt_bdint_path_solve_sam_singloc(path, solseries, sing_start); % Gets the singularity location
%[~, singcoeff(3,k)] = qt_bdint_path_solve_sam_singcoeff(path, solseries, sing_start, -singloc(k)); % Gets the path and the singularity coefficient
singang(3,k) = atan(imag(singloc(3,k))/real(singloc(3,k)));
end


%======== RIGHT DIAGONAL ========%

path = get_path('right_diag_singularity');
sing_start = 12;

for k=1:M
% Initial guess - last two are mu and B (in that order!)
x0 = zeros(1,N-1);
x0(N) = 1;
x0(N+1) = 0.1;
ep = eps(k);
est = zeros(1,N+1);

% Solve the ODE
mytol = 1e-13;
options = optimset('Display', 'iter', 'TolFun', mytol, 'TolX', mytol);
fwd = @(x)series_bern(x, ep, N, kap, grav_coeff);
[est,fval] = fsolve(fwd,x0,options); % Call solver

% Save all the data
solseries = [];
solseries.N = N;
solseries.an = est(1:N-1);
solseries.mu = est(N);
solseries.B = est(N+1);
solseries.ep = ep;
solseries.kap = kap;
solseries.ycrest = sum(est(1:N-1)./(2*pi*(1:N-1)));

[~, singloc(4,k)] = qt_bdint_path_solve_sam_singloc(path, solseries, sing_start); % Gets the singularity location
%[~, singcoeff(4,k)] = qt_bdint_path_solve_sam_singcoeff(path, solseries, sing_start, -singloc(k)); % Gets the path and the singularity coefficient
singloc(4,k) = -singloc(4,k)

end

%{
%======== LEFT DIAGONAL (ABOVE SHEET) ========%

path = get_path('left_diag_singularity2');
sing_start = 20;

for k=1:M
% Initial guess - last two are mu and B (in that order!)
x0 = zeros(1,N-1);
x0(N) = 1;
x0(N+1) = 0.1;
ep = eps(k);
est = zeros(1,N+1);

% Solve the ODE
mytol = 1e-13;
options = optimset('Display', 'iter', 'TolFun', mytol, 'TolX', mytol);
fwd = @(x)series_bern(x, ep, N, kap, grav_coeff);
[est,fval] = fsolve(fwd,x0,options); % Call solver

% Save all the data
solseries = [];
solseries.N = N;
solseries.an = est(1:N-1);
solseries.mu = est(N);
solseries.B = est(N+1);
solseries.ep = ep;
solseries.kap = kap;
solseries.ycrest = sum(est(1:N-1)./(2*pi*(1:N-1)));

[~, singloc(5,k)] = qt_bdint_path_solve_sam_singloc(path, solseries, sing_start); % Gets the singularity location
%[~, singcoeff(5,k)] = qt_bdint_path_solve_sam_singcoeff(path, solseries, sing_start, -singloc(k)); % Gets the path and the singularity coefficient

end

%======== LEFT DIAGONAL (ABOVE 2 SHEETs) ========%

path = get_path('left_diag_singularity3');
sing_start = 28;

for k=1:M
% Initial guess - last two are mu and B (in that order!)
x0 = zeros(1,N-1);
x0(N) = 1;
x0(N+1) = 0.1;
ep = eps(k);
est = zeros(1,N+1);

% Solve the ODE
mytol = 1e-13;
options = optimset('Display', 'iter', 'TolFun', mytol, 'TolX', mytol);
fwd = @(x)series_bern(x, ep, N, kap, grav_coeff);
[est,fval] = fsolve(fwd,x0,options); % Call solver

% Save all the data
solseries = [];
solseries.N = N;
solseries.an = est(1:N-1);
solseries.mu = est(N);
solseries.B = est(N+1);
solseries.ep = ep;
solseries.kap = kap;
solseries.ycrest = sum(est(1:N-1)./(2*pi*(1:N-1)));

[~, singloc(6,k)] = qt_bdint_path_solve_sam_singloc(path, solseries, sing_start); % Gets the singularity location
%[~, singcoeff(6,k)] = qt_bdint_path_solve_sam_singcoeff(path, solseries, sing_start, -singloc(k)); % Gets the path and the singularity coefficient

end

%}

%======== EXTRA UHP LEFT SINGULARITY (NEEDS BETTER NAME) ========%

path = get_path('extra_singularity_left');
sing_start = 12;

for k=1:M
% Initial guess - last two are mu and B (in that order!)
x0 = zeros(1,N-1);
x0(N) = 1;
x0(N+1) = 0.1;
ep = eps(k);
est = zeros(1,N+1);

% Solve the ODE
mytol = 1e-13;
options = optimset('Display', 'iter', 'TolFun', mytol, 'TolX', mytol);
fwd = @(x)series_bern(x, ep, N, kap, grav_coeff);
[est,fval] = fsolve(fwd,x0,options); % Call solver

% Save all the data
solseries = [];
solseries.N = N;
solseries.an = est(1:N-1);
solseries.mu = est(N);
solseries.B = est(N+1);
solseries.ep = ep;
solseries.kap = kap;
solseries.ycrest = sum(est(1:N-1)./(2*pi*(1:N-1)));

[~, singloc(7,k)] = qt_bdint_path_solve_sam_singloc(path, solseries, sing_start); % Gets the singularity location
%[~, singcoeff(6,k)] = qt_bdint_path_solve_sam_singcoeff(path, solseries, sing_start, -singloc(k)); % Gets the path and the singularity coefficient
singloc(7,k) = -singloc(7,k);
end

%======== EXTRA UHP LEFT SINGULARITY (NEEDS BETTER NAME) ========%

path = get_path('extra_singularity_right');
sing_start = 12;

for k=1:M
% Initial guess - last two are mu and B (in that order!)
x0 = zeros(1,N-1);
x0(N) = 1;
x0(N+1) = 0.1;
ep = eps(k);
est = zeros(1,N+1);

% Solve the ODE
mytol = 1e-13;
options = optimset('Display', 'iter', 'TolFun', mytol, 'TolX', mytol);
fwd = @(x)series_bern(x, ep, N, kap, grav_coeff);
[est,fval] = fsolve(fwd,x0,options); % Call solver

% Save all the data
solseries = [];
solseries.N = N;
solseries.an = est(1:N-1);
solseries.mu = est(N);
solseries.B = est(N+1);
solseries.ep = ep;
solseries.kap = kap;
solseries.ycrest = sum(est(1:N-1)./(2*pi*(1:N-1)));

[~, singloc(8,k)] = qt_bdint_path_solve_sam_singloc(path, solseries, sing_start); % Gets the singularity location
%[~, singcoeff(6,k)] = qt_bdint_path_solve_sam_singcoeff(path, solseries, sing_start, -singloc(k)); % Gets the path and the singularity coefficient

end


%% Plotting stuff

figure(3);
plot(eps, central_diff(log(abs(singloc(1,:))), eps));
hold on
plot(eps, central_diff(log(abs(singloc(3,:))), eps));
title('gradient of log plot for singularities');
legend('First singularity', 'left diagonal singularity');

figure(1);
%label('ep', 'abs(z)');
plot(eps, abs(singloc(1,:)))
hold on
plot(eps, abs(singloc(2,:)))
plot(eps, abs(singloc(3,:)))
plot(eps, abs(singloc(4,:)))
plot(eps, abs(singloc(5,:)))
plot(eps, abs(singloc(6,:)))
plot(eps, abs(singloc(7,:)))
legend('First singularity', 'unphysical fluid singularity 1', 'left diagonal singularity (LHP)', 'right diagonal singularity (LHP)' ,'left diag. sing. (sheet above)', 'left diag. up again', 'UHP left diagonal');


figure(2);
plot(real(singloc(1,:)), imag(singloc(1,:)));
hold on
plot(real(singloc(2,:)), imag(singloc(2,:)));
plot(real(singloc(3,:)), imag(singloc(3,:)));
plot(real(singloc(4,:)), imag(singloc(4,:)));
plot(real(singloc(5,:)), imag(singloc(5,:)));
plot(real(singloc(6,:)), imag(singloc(6,:)));
plot(real(singloc(7,:)), imag(singloc(7,:)));
plot(real(singloc(8,:)), imag(singloc(8,:)));
legend('First singularity', 'unphysical fluid singularity 1', 'left diagonal singularity', 'right diagonal singularity' ,'left diag. sing. (sheet above)', 'left diag. up again');


figure(3);
plot(eps, singang(3,:));
legend('First singularity', 'unphysical fluid singularity 1', 'left diagonal singularity', 'right diagonal singularity' ,'left diag. sing. (sheet above)', 'left diag. up again');

figure(4);
plot(real(singloc(1,:)), imag(singloc(1,:)));
hold on
plot(real(singloc(2,:)), imag(singloc(2,:)));
plot(real(singloc(3,:)), imag(singloc(3,:)));
plot(real(singloc(4,:)), imag(singloc(4,:)));
plot(real(singloc(7,:)), imag(singloc(7,:)));
plot(real(singloc(8,:)), imag(singloc(8,:)));
legend('First singularity', 'unphysical fluid singularity', 'left LHP diagonal singularity', 'right LHP diagonal singularity', 'left UHP diagonal singularity', 'right UHP diagonal singularity');


