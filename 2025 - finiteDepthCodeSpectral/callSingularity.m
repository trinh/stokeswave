clear


%% Script to get singularity for diffrent wave energies and mass.

%Generate leading order real wave.

L=0.5;%Half the x period.
N=1*1024;%Number of gridpoints (or Fourier coefficients).
dxi=2*L/N;
xi=-L:dxi:L-dxi;
n=[0:N/2-1 0 -N/2+1:-1];%Fourier summation index in FFT
Mk=1i*n*pi/L;%Fourier multiplier for differentiation


%% Generate leading order solution
%data = load("testN=150.mat");
%MV = linspace(0.08,0.1,5); %Mass
%MV = data.MMat;
%AmpV = linspace(3,0.1,5)/10000;
%AmpV = data.EMat;

%% Guve energy and mass of wave.
%This can be given as a vector if so we will loop over all mass and
%energies

MV = 0.4628;
AmpV = 5.5507e-4;
%% Loop over doffrent energies and masses
Yout = [];
fsing = [];
hvec = [];
Omegavec = [];
for count2 = 1
    
    for count1 = 1:length(AmpV)
        %%Get real soloutions
        Amp = AmpV(count1);
        M = MV(count1);
        E = Amp;
        B = 0; %Bond number
        C = M;
        Yint = M+Amp*cos(2*pi*xi);
        Y = Yint;
        Yhat = fft(Y);
        Mt = 1i*coth(2*pi*mean(Y)*n);%Fourier multiplier for the T integral.
        X = xi-real(ifft(Mt.*Yhat));
        F = sqrt((1+(2*pi)^2*B)/((2*pi)*coth(2*pi*M)));
        
        %Converting in to tau theta varibles
        X1 = 1-real(ifft(Mt.*Mk.*Yhat));
        Y1 = real(ifft(Mk.*Yhat));
        Omegavec = [Omegavec; log(1./(X1+1i*Y1))];

        %% Newton Iteration. %%
        %iAmp = 1;%Iterate on amplitude.
        iAmp = 2;%Iterate on energy.
        iEig = 0;%Froude number eigenvalue
        % iEig = 1;%Bond number eigenvalue
        [Yreal,C,F,B,M,err,fval] = NewtonSolver_Bern_Fsolve(Y,C,F,B,M,E,iAmp,iEig);
        %----------------------%
        h = mean(Yreal);
        Yout = [Yout;Yreal];
        
        %%
        
        % % Linear solution. %
        %--- Key parameters ---%
        if B ~=0 %Ensure input data is linear solotuion
            disp('B must be equal to zero')
            return
        end
        Amp = E;
        %----------------------%
       
        Yhat = fft(Yreal);
        Mt = 1i*coth(2*pi*mean(Yreal)*n);%Fourier multiplier for the T integral.
        %% Set up the Intergartion mesh.
        
        %Set up mesh
        Nxi = N; Neta = 100;
        xi_mesh = linspace(-0.5, 0.5-dxi, Nxi);
        eta_mesh = linspace(0, 2*h, Neta); %Code is only valide for 0<2h
        
        % Make mesh
        [XI, ETA] = meshgrid(xi_mesh, eta_mesh);
        ZETA = XI + 1i*ETA;
        % Initialise ODE
        Ymat = 0*ZETA;
    
        mytol = 1e-10;
        Ifunc = @(Z) getI(Yreal,Z,N);
        myfunc = @(Z,Y) getDY(Y,F,C,Z,Ifunc);
        sfun = @(t) eta_mesh(1) + (eta_mesh(end) - eta_mesh(1))*t; %We need sfunc for correc stop conditons 
        opts = odeset('Events', @(t,y) steepEvent_time(t,y,myfunc,sfun),'RelTol',mytol, 'AbsTol', mytol, 'Stats','off');
        %% Get Singularity intergrate up imaginary axis until we hit singularity
      
        path = 0+[eta_mesh(1), eta_mesh(end)]*1i; 

        Yinit = interp1(xi, Yreal, 0);

        [Svec, out] = get_ode_line(path, Yinit, myfunc, opts);
        %Yvec = interp1(imag(Svec), out, eta_mesh);
        fsing = [fsing; imag(Svec(end))];
        hvec = [hvec, h];
        %% Ploting stuff
        %plot
        % figure(10)
        % hold on
        % plot(xi,Yreal,Color='red')
        % figure(11)
        % hold on
        % semilogy([0:N/2-2],abs(fft(Yreal(1:N/2-1))),'blue')
        % xlim([0 N/2-2])
        % ylim([10^(-17) 1000])
    end
end
%%
% figure(1)
% plot(data.r0Vec,data.fstarMat,'o-')
% hold on 
% plot(data.r0Vec(2:end),fsing,'o-')
% plot(data.r0Vec(2:end),2*hvec,'o-')
% ylim([0.05 0.2])
% 
% figure(2)
% 
% plot(xi,Yout,'blue')
% hold on
% plot(data.xi,data.Yy,'red')
% 
% figure(3)
% hold on
% plot(100*ones(15,1),abs(data.Yy(2:end,512)-Yout(:,N/2)),'o')



function [value, isterminal, direction] = steepEvent_time(t, y, myfunc,sfun)
    %Stop condtion for ode.
    %Stop when the dardient of dy/dx is to steep.
    z = sfun(t);
    dydt = myfunc(z, y);
   
    value = 1/abs(dydt)-5e-7;

    isterminal = 1;                 % stop
    direction = -1;                 % only when steepness increases past threshold
end
