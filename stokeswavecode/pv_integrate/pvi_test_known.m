% We perform some tests to find a scheme that can integrate the PVI
% We are interested in seeing how the integration proceeds as we approach
% the axis
%
% Note this may not work particularly well near the edges of the
% integration (-0.5, 0.5) in our Stokes wave

function pvi_test_known

% Test singularities
w = -0.4 + 1e0i; % 0.0000302159 + 0.0000219531 I
w = -0.4 + 1e-1i; % 0.00863202 + 0.00627153 I
w = -0.4 + 1e-2i; % 0.015195 + 0.0110398 I
w = -0.4 + 1e-3i; % 0.016079 + 0.0116821 I
w = -0.4 + 1e-5i; % 0.016079 + 0.0116821 I
w = -0.4 + 1e-11i; % 0.0161803
w = -0.4 + 0i; % 0.0161803

% Approximate Stokes wave for ep = 0.05
a1 = 0.02;
f = @(u) imag(1 + a1*exp(-2i*pi*u));
F = @(u) f(u).*cot(pi*(u - w));
disp('========================================');
residue = f(real(w));
disp(['Residue on the axis gives = ', num2str(residue)]);
disp('========================================');

N = 1000;
phi = linspace(-0.5, 0.5, N);

tic;
I = trapz(phi, F(phi));
t = toc;
disp(['Trapezoid = ', num2str(I), ' (time = ', num2str(t), ')']);

tic;
I = simps(phi, F(phi));
t = toc;
disp(['Simpsons = ', num2str(I), ' (time = ', num2str(t), ')']);

tic
I = integral(F, -0.5, 0.5);
t = toc;
disp(['Integral = ', num2str(I), ' (time = ', num2str(t), ')']);

tic
[I, collapse] = pvi_gauss(F, w, -0.5, 0.5);
% If things were collapsed, we'll need to add in the residue manually. Note
% it depends on whether we approach from top or bottom
if collapse == 1    
    I = I + 1i*residue;
end
t = toc;
disp(['PVI-Gauss = ', num2str(I), ' (time = ', num2str(t), ')']);
