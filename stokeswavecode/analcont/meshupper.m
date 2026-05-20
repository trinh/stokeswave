clear
close all

% =====================================
%         Load in the solution
F = load('../../stokeswavedata/series_N100_ep0.2.mat');
solseries = F.solseries;
N = solseries.N;
mu = solseries.mu;
B = solseries.B;
an = solseries.an;
% =====================================

% Initial mesh for UHP
M = 100; % No. of meshpoints
phi = linspace(-0.5, 0.5, M);
psiUHP = linspace(0, 1, M);
[Phi,PsiUHP] = meshgrid(phi,psiUHP);
fUHP = Phi+1i*PsiUHP;

% Form the initial conditions
numbranch = 2;
gaxis = zeros(numbranch, length(phi));
gaxis(1, :) = zfree(phi, solseries) + zfree(-phi, solseries);

% Initial loop
psitop = 2i;
path = [phi(end), ...
    phi(end) + psitop, ...
    phi(1) + psitop, ...
    phi(1)];

[ftmp, gtmp] = gbern_contour(path, solseries);

ga = 0*phi;
ga(1) = gtmp(end);
for j = 1:length(phi)-1
    [~, tmp2] = gbern_contour([phi(j), phi(j+1)], solseries, ga(j));
    ga(j+1) = tmp2(end);
end

gaxis(2,:) = ga;

meshdat = [];
for k = 1:numbranch
    
    gUHP = zeros(M,M);
    
    for j=1:length(phi)
        gam = @(s) phi(j) + 1i*s;
        dgam = 1i;
        
        fwd = @(s,g)gbern_lower2upper(gam(s), g, solseries)*dgam;
        
        %initial condition
%         g0 = zfree(phi(j), solseries) + zfree(-phi(j), solseries);
        
        [a, gUHP(:,j)] = ode45(fwd,psiUHP, gaxis(k, j));
    end
    
    %Now g(f) = z(f) + z(-f) so we need to subtract z(-f) = zminusf
    zminusf = zfree(-fUHP, solseries);
    
    % Therefore z(f) (where f is in UHP)
    zUHP = gUHP - zminusf;
    
    pUHP = ones(size(fUHP));
    for j = 1:N-1
        pUHP = pUHP + an(j).*exp(2*pi*1i*j*fUHP);
    end
    gpUHP = mu./(2*pUHP.*(mu*B + 1i*pi*gUHP)) - pUHP;
    zpUHP = gpUHP + pUHP;
    
    % Save the mesh
    
    meshdat.solseries = solseries;
    meshdat.UHP(k).f = fUHP;
    meshdat.UHP(k).g = gUHP;
    meshdat.UHP(k).gp = gpUHP;
    meshdat.UHP(k).p = pUHP;
    meshdat.UHP(k).z = zUHP;
    meshdat.UHP(k).zp = zpUHP;
    
    
end

%plot the in fluid solution to match up
psiLHP = linspace(-1,0,M);
[~,PsiLHP] = meshgrid(phi,psiLHP);
fLHP = Phi+1i*PsiLHP;
zLHP = zfree(fLHP, solseries);
meshdat.LHP.f = fLHP;
meshdat.LHP.z = zLHP;

% name = ['meshdat_M', num2str(M), '_ep', num2str(solseries.ep), '.mat'];
% save(['../data/', name], 'meshdat');
