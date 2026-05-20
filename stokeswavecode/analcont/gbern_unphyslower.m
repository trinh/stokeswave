clear
close all

% =====================================
%         Load in the solution
F = load('data/series_N200_ep0.05.mat');
solseries = F.solseries;

F = load('data/meshdat_M100_ep0.05.mat');
meshdat = F.meshdat;
% =====================================

fLHP = meshdat.LHP.f;
phi = real(fLHP(1,:));
psi = imag(fLHP(:,1));
gLHP2 = 0*fLHP;
branchselect = 1;
gselect = 2;
pU = meshdat.UHP(branchselect).zp;
fU = meshdat.UHP(branchselect).f;

ind = phi < 0;
fUHP.left = fU(:,ind);
pUHP.left = pU(:, ind);
ind = phi >= 0;
fUHP.right = fU(:, ind);
pUHP.right = pU(:, ind);

% Initial condition on g
gaxis = meshdat.UHP(gselect).g(1,:);


for j = 1:length(phi)
    fwd = @(f, g) 1i*gbern_upper2lower(phi(j) + 1i*f, g, solseries, fUHP, pUHP);
    
    %initial condition - NEED TO TAKE UPPER FREE SURFACE
    g0 = gaxis(j);
    [temp1, temp2] = ode45(fwd, flipud(psi), g0); 
    gLHP2(:,j) = flipud(temp2);    
%     if j == 45 || j == 49 || j == 50 || j == 51
%         keyboard
%     end
end

zLHP2 = gLHP2 - rot90(meshdat.UHP(branchselect).z, 2);

%%
comp = 5; % Compression

figure(1); figshift; %clf(1); 
hold on
plotme(fLHP, imag(gLHP2), 'g', comp);
plotme(fU, imag(meshdat.UHP(1).g), 'b', comp);
plotme(fU, imag(meshdat.UHP(2).g), 'r', comp);

figure(2); figshift; %clf(2);
hold on
plotme(fLHP, real(zLHP2), 'g', comp);
plotme(fU, real(meshdat.UHP(1).z), 'b', comp);
plotme(fU, real(meshdat.UHP(2).z), 'r', comp);

figure(3); figshift; %clf(2);
hold on
plotme(fLHP, imag(zLHP2), 'g', comp);
plotme(fU, imag(meshdat.UHP(1).z), 'b', comp);
plotme(fU, imag(meshdat.UHP(2).z), 'r', comp);

