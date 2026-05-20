% Code to generate analytic continuation paths along arbitrary contours 

clear
close all

% =====================================
%         Load in the solution
r0 = 0.7;
ep = sqrt(0.4);
wave = VBStokesWave(r0, ep);
Q = -log(r0)/(2*pi); %For easy postmortum access
% =====================================


% =====================================
%            Path input
% path = get_path('round.0.A');
P = get_path('box1pass');
path.singstart = P.singstart;
path.f = split_path(P.f,r0);
[D,L,UD]=region_indicator(path.f,r0);
path.L=L;
path.UD=UD;

pathdat = qt_bdint_path_solve(path, wave);
% pathdat = qt_bdint_path_solve(path, solseries, 'singcoef', [1, 0.17555i]);


%%
omega = pathdat.tau - 1i*pathdat.theta;
W = exp(omega);
[phi,psi,FP]=SeriesData(wave, pathdat);

% omega2 = pathdat.tau2 - 1i*pathdat.theta2;
% W2 = [NaN*ones(length(W)), exp(omega2)];
% 
% omega3 = pathdat.tau3 - 1i*pathdat.theta3;
% W3 = [NaN*ones(2*length(w)), exp(omega3)];


figure(5); figshift; hold all
% plot3(real(pathdat.f),imag(pathdat.f), imag(pathdat.omega), '-', 'LineWidth', 2);
% plot3(real(pathdat.f),imag(pathdat.f), imag(W), '-', 'LineWidth', 2);
plot3(real(pathdat.f),imag(pathdat.f), imag(W), '-', 'LineWidth', 2);
xlabel('\phi'); ylabel('\psi');
hold on
surf(phi,psi,imag(FP),'linestyle', 'none');
alpha 0.5




% Run the following code to plot
% handles = qt_bdint_path_plotter(path, pathdat, 'PatchWidth', 0.02, 'NumArc', 150); 


figure(6); figshift; hold all
% plot3(real(pathdat.f),imag(pathdat.f), imag(pathdat.omega), '-', 'LineWidth', 2);
% plot3(real(pathdat.f),imag(pathdat.f), imag(W), '-', 'LineWidth', 2);
plot3(real(pathdat.f),imag(pathdat.f), real(W), '-', 'LineWidth', 2);
xlabel('\phi'); ylabel('\psi');
surf(phi,psi,real(FP),'linestyle', 'none');
alpha 0.5