close all
clear all

N=40;
ep = 0.8;
kap = 0;
grav_coeff = 1;

x0 = zeros(1,N-1);
x0(N) = 1;
x0(N+1) = 0.1;

% Solve the ODE
mytol = 1e-13;
options = optimset('Display', 'iter', 'TolFun', mytol, 'TolX', mytol);
fwd = @(x)series_bern(x, ep, N, kap, grav_coeff);
[est,fval] = fsolve(fwd,x0,options); % Call solver

%%
% Wave height
z_crest = 0;
z_trough = 0.5;

for k=1:N-1
    z_crest = z_crest + ((1i*est(k))/(2*pi*k));
    z_trough = z_trough + ((1i*est(k))/(2*pi*k))*(exp(-2*k*pi*1i*0.5));
end

height = imag(z_crest-z_trough)


%%
k2 = linspace(1,(N-1),N-1);
test_alpha = ((-1)./(2.*(k2))).*log(k2)-((log(abs(est(k2))))./(k2));
test_alpha2 = (1./((k2))).*((-0.5).*log((k2)) - log(abs(est(k2))));
plot(k2,test_alpha);
hold on
plot(k2,test_alpha2);

%{
% Compare the asymptotic behaviour
x = linspace(1,N-1,N-1);
asym = ((x).^(-1.5)).*exp(-k*0.01);
figure(1);
hold on
plot(x, est(1:N-1));
plot(x, asym);
%}

spectral_estimate = (1/((N-1)))*((-0.5)*log((N-1)) - log(abs(est(N-1))))
radius_estimate = (1/(2*pi))*log(1/((abs(est(N-1)))^(1/(N-1))))