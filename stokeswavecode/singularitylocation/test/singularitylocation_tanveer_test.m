clear
close all


% Test with known function: G(f) = (f-f0)^1/2
% Gp = (1/2)G^(-1)

f0 = 0.45*1i; % Singularity location


% Contour containing singularity: go around twice for 1/2 branch point.
path = [0.5, 0.5+1i, -0.5+1i, -0.5, 0.5, 0.5+1i, -0.5+1i, -0.5, 0.5];

[f, G, integral] = odecontour(path,f0);

figure(1);
plot3(real(f),imag(f),imag(G));


(integral)*(1/(2*pi*1i))


%{
location = (1/(2*pi*1i))*integral
height = pi*(imag(zfree(0,solseries) - zfree(0.5,solseries)))
%}




