% 2026 edit: 
%  - Testing with Henry to check if boundary integral returns same solution 
%  - Added flag to not overwrite


clear
close all

N = 51; % Should be odd
ep = 0.7746; % Amplitude
r0 = 0.3; % Depth parameter

% Current analytic continuaton test
% r0 = 0.5;
% ep = sqrt(0.1);

saveme = 0;

name = ['bdint','_N', num2str(N),'_r', num2str(r0),'_ep', num2str(ep),'.mat'];

if exist('name', 'file') == 0
    disp('No existing data...solving new...')
  %if no solution data exists, solve as usual 
    phis =  linspace(-0.5,0.5, N);
    wave = VBStokesWave(r0, ep);
    z = [];
    zp = [];
        for j=1:length(phis)
            [Z, Zp] = getZValues(wave, phis(j));
            zp(j)= Zp;
            z(j) = Z;
        end
        
    bdint_data = [];
    bdint_data.zp = zp;
    delta_guess = real(zp);
    beta_guess = imag(zp);
    
    x0 = zeros(1,2*N-2);
    x0(2*N-1) = wave.B;
    x0(2*N) = wave.mu;

    delta_guess_interp = interp1(linspace(-0.5,0.5,N), delta_guess, linspace(-0.5,0.5,N));
    beta_guess_interp = interp1(linspace(-0.5,0.5,N), beta_guess, linspace(-0.5,0.5,N));
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
    %integrate dz/dphi using cumtrapz to get z
    phi2 = linspace(-0.5,0.5,N);
    %boundary integral result
    x = cumtrapz(phi2,delta);
    y = cumtrapz(phi2,beta);
    %series truncation result
    x2 = cumtrapz(phi2,delta_guess_interp); 
    y2 = cumtrapz(phi2,beta_guess_interp);
    
    %%
    bdint_data.x = x;
    bdint_data.y = y;
    bdint_data.x2 = x2;
    bdint_data.y2 = y2;
    name = ['./../bdint_data/', 'bdint','_N', num2str(N),'_r', num2str(r0),'_ep', num2str(ep),'.mat'];

    if saveme == 1
        save(name, 'bdint_data')
    end

else
    %load solution data if it exists
    name = ['../bdint_data/', 'bdint','_N', num2str(N),'_r', num2str(r0),'_ep', num2str(ep),'.mat'];
    L = load(name);
    loaded_data = L.bdint_data;
    zp = loaded_data.zp;
    
 
    x = loaded_data.x;
    y = loaded_data.y;
    x2 = loaded_data.x2;
    y2 = loaded_data.y2;
end
 
%plot the surface
%%
figure(2);
plot(x,y, 'k'); hold all
plot(x2,y2, 'ko');
axis equal

legend('Boundary integral', 'Series truncation');

