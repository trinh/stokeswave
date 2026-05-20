clear
close all

N = 101; % Should be odd
ep = 0.8; % Amplitude
r0 = 0; % Depth parameter

% Load initial guess from infinite depth
load('point3initialguessN51');

x0 = zeros(1,2*N-2);
x0(2*N-1) = 1;
x0(2*N) = 1;

delta_guess_interp = interp1(linspace(-0.5,0.5,51), delta_guess, linspace(-0.5,0.5,N));
beta_guess_interp = interp1(linspace(-0.5,0.5,51), beta_guess, linspace(-0.5,0.5,N));
x0(1:N) = delta_guess_interp;
x0(N+1:2*N-2) = beta_guess_interp(2:N-1);

tic
options = optimset('Display', 'iter');
fwd = @(X)boundaryintegral(X, N, ep, r0);
[sol,fval] = fsolve(fwd,x0,options); % Call solver
toc

delta = sol(1:N);
beta(1) = 0;
beta(N) = 0;
for k=2:N-1
    beta(k)=sol(N+k-1); 
end
B = sol(2*N-1);
mu = sol(2*N);

phi2 = linspace(-0.5,0.5,N);
x = cumtrapz(phi2,delta);
y = cumtrapz(phi2,beta);
x2 = cumtrapz(phi2,delta_guess_interp);
y2 = cumtrapz(phi2,beta_guess_interp);

%%

figure(1);
plot(x,y, '-'); hold all
plot(x2,y2, 'o');
axis equal

legend('Boundary integral', 'Series truncation');


