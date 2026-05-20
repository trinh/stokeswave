% Analytically continue using the boundary integral formulation

clear
close all

% =====================================
%         Load in the solution
F = load('../data/bdint_N201_ep0.05.mat');
solbdint = F.solbdint;
ep = solbdint.ep;
F = load('../data/series_N100_ep0.05.mat');
solseries = F.solseries;
% =====================================

% So we introduce an 'a' flag and also a inLHR flag. 
% 'a' marks the initial direction of analytic continuation (+1 for UHP)
% inLHR is a switch that flips on when we re-enter the LHR. 
% It is unclear whether both are the same
a = -1;
inLHP = 0;
going = 'down';

psimax = 2;
psimin = 1e-3;

phi = solbdint.phi;

M = length(phi);
switch going
    case 'up'
        psi = linspace(psimin, psimax, M);
    case 'down'
        psi = linspace(-psimax, -psimin, M);
    otherwise 
        error('Check sign of dir');
end
[PHI, PSI] = meshgrid(phi, psi);
wmat = PHI + 1i*PSI;

ymat = 0*wmat; ypmat = 0*wmat;
xmat = 0*wmat; xpmat = 0*wmat;

for j = 1:length(phi)
    fwd = @(psi,Y) bdfunc_complex(phi(j) + 1i*psi, Y, solbdint, solseries, a, inLHP);
    % Note d/dpsi = i d/dpsi
    fwdpath = @(psi, Y) 1i*fwd(psi, Y);
    switch going
        case 'up'
            [T, Y0] = ode45(fwdpath, psi, [solbdint.x(j), solbdint.y(j)]);   
        case 'down'
            [T, Y0] = ode45(fwdpath, fliplr(psi), [solbdint.x(j), solbdint.y(j)]);   
            T = flipud(T);
            Y0 = flipud(Y0);
        otherwise 
            error('Check sign of dir');
    end
    xmat(:,j) = Y0(:,1);
    ymat(:,j) = Y0(:,2);    
    
    
    % Send it back in to get the derivatives (note annoying shift of sizes)
    % I don't know why this does not work well
%     keyboard
    Y0p = fwd(T.', Y0.');
    xpmat(:, j) = Y0p(1,:);
    ypmat(:, j) = Y0p(2,:);    
%     keyboard
%     xpmat(:,j) = -1i*central_diff(xmat(:,j), psi);
%     ypmat(:,j) = -1i*central_diff(ymat(:,j), psi);    
end

% For some reason, the first element always gets screwed up
% We extend the mesh to include psi = 0, and include the initial conditions
switch going
    case 'up'
        keyboard
    case 'down'
        % The first element always gets screwed up!
%         xpmat(end,:) = central_diff(xmat(end,:), phi);
%         ypmat(end,:) = central_diff(ymat(end,:), phi);
        
        wmat(end+1,:) = phi;
        xmat(end+1, :) = solbdint.x;
        ymat(end+1, :) = solbdint.y;
        xpmat(end+1, :) = central_diff(solbdint.x, phi);
        ypmat(end+1, :) = central_diff(solbdint.y, phi);
%         xpmat(end+1, :) = central_diff(xmat(end,:), phi);
%         ypmat(end+1, :) = central_diff(ymat(end,:), phi);
end
zmat = xmat + 1i*ymat;
zpmat = xpmat + 1i*ypmat;

meshdatbd = [];
meshdatbd.wmat = wmat;
meshdatbd.xmat = xmat;
meshdatbd.xpmat = xpmat;
meshdatbd.ymat = ymat;
meshdatbd.ypmat = ypmat;
meshdatbd.B = solbdint.B;
meshdatbd.mu = solbdint.mu;
meshdatbd.ep = solbdint.ep;
meshdatbd.a = a;
meshdatbd.inLHP = inLHP;
meshdatbd.going = going;

%%
name = ['./../data/meshdatbd_ep', num2str(ep), '.mat'];
save(name, 'meshdatbd');

%%
% figure(1); clf(1); figshift
% subplot(1, 2, 1)
% surf(real(w), imag(w), real(xmat));
% 
% 
% subplot(1, 2, 2)
% surf(real(w), imag(w), imag(xmat));
% 
% 
% 
% figure(2); clf(2); figshift
% subplot(1, 2, 1)
% surf(real(w), imag(w), real(ymat));
% 
% 
% subplot(1, 2, 2)
% surf(real(w), imag(w), imag(ymat));

figure(3); clf(3); figshift
% subplot(1, 2, 1)
surf(real(wmat), imag(wmat), real(xmat));
shading interp

% subplot(1, 2, 2)
% surf(real(w), imag(w), imag(zmat));



%{
rowtest = length(psi);
ztop = zmat(rowtest,:);
zptop = central_diff(ztop, phi);
bernerror = 1./(2*abs(zptop).^2) + 2*pi/solbdint.mu.*ymat(end,:) - solbdint.B;
figure(3); clf(3);
plot(phi, abs(bernerror));
%}