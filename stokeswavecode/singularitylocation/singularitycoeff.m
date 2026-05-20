% =====================================
%         Load in the solution
F = load('../../stokeswavedata/series_N100_ep0.3.mat');
solseries = F.solseries;
% =====================================

p=0.4;
M = 25; % For the interpolation
path = [p, p+1i, -p+1i, -p,p, p+1i, -p+1i, -p,p];

[f2,g2,gp2,zp2] = gbern_contour(path,solseries);

integral = 0;
figure(1);
hold on

for k=1:length(path)-1
    dgam = (1i)^(k);
    path2 = path(k:k+1);
    if k == 1
        [f,g,gp,zp] = gbern_contour(path2,solseries);
        z = g-zfree(-f,solseries);
    else
        g0 = g(end);
        [f,g,gp,zp] = gbern_contour(path2,solseries,g0);
        z = g-zfree(-f,solseries);
    end
    
    if mod(k,4) == 0
        t = real(f);
    elseif mod(k,4) == 1
        t = imag(f);
    elseif mod(k,4) == 2
        t = -real(f);
    else 
        t = -imag(f);
    end
    
    integrand = z.*(f-0.5*1i);
    integrand = integrand*dgam;
    integral = integral + trapz(t, integrand);
    plot3(real(f),imag(f),imag(z));
end

coefficient = (1/(2*pi*1i))*integral



