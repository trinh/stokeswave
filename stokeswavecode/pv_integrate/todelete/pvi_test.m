% We perform some tests to find a scheme that can integrate the PVI

function pvi_test

% =====================================
%         Load in the solution
F = load('./../data/series_N100_ep0.05.mat');
solseries = F.solseries;
% =====================================

mytol = 1e-8;

N = 1000;
phi = linspace(-0.5, 0.5, N);

w = -0.4 + 1e-8i;
w = -0.4 + 1e-2i;
w = -0.4 + 0*1e-20i;

I = trapz(phi, fullfunc(phi, w));
disp(['Trapezoid = ', num2str(I)]);

I = simps(phi, fullfunc(phi, w));
disp(['Simpsons = ', num2str(I)]);

tic
I = integral(@(f)fullfunc(f, w), -0.5, 0.5);
disp(['Integral = ', num2str(I)]);
toc

% Try the PVI routine
% Need to provide derivatives
f = @(u) imag(zpfree(u, solseries));
fp = @(u) imag(zppfree(u, solseries));
h = @(u) tan(pi*(u - w));
hp = @(u) pi*sec(pi*(u - w)).^2;
hpp = @(u) 2*pi^2*sec(pi*(u - w)).^2.*tan(pi*(u - w));
tic
I = pvi_nyiri(f, fp, h, hp, hpp, w, -0.5, 0.5);
disp(['Integral Nyiri = ', num2str(I)]);
toc

    function g = fullfunc(f, w0)
        [~, zp] = zfree(f, solseries);
        g = imag(zp).*cot(pi*(f - w0));        
    end
end