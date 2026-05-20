% Non-monotonic relationship between ep and mu
clear all 
close all

N = 50; 
M = 100; %ep iterations

eps = linspace(0.6,0.95,M);
est = zeros(N+1,M);

for k=1:M
ep = eps(k);
kap = 0;
grav_coeff = 1;

% Initial guess - last two are mu and B (in that order!)
if k == 1
    x0 = zeros(1,N-1);
    x0(N) = 1;
    x0(N+1) = 0.1;
else
    x0 = est(:,k-1);
    x0 = transpose(x0);
end

% Solve the ODE
mytol = 1e-13;
options = optimset('Display', 'iter', 'TolFun', mytol, 'TolX', mytol);
fwd = @(x)series_bern(x, ep, N, kap, grav_coeff);
[est(:,k),fval] = fsolve(fwd,x0,options); % Call solver

end

%%
% Plot ep against mu
 
figure();
plot(eps,est(N,:));
