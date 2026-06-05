%Code manually serch for singularity.
%We can give a custome path to investigate the sheet structor

clear


%We requie a data file which contains
%M  'mass'
%F  'Froude number'
%Yreal 'Real solotuion for the free surface'
%Y0     'Solotion of Y on real axis for non primary sheets can be gotten
%callManualPath

%laod data
data = load("A0Deep.mat");

%Set up paramters and calculate X,Xxi,Yxi,h
%Todo this can be put as a sperate function, and data structor. 
L=0.5;%Half the x period.
N= data.N;

dxi=2*L/N;
xi=-L:dxi:L-dxi;
n=[0:N/2-1 0 -N/2+1:-1]; %Fourier summation index in FFT
Mk=1i*n*pi/L;

M = data.M; 
F = data.F;
C = data.M;     % Bernoulli constant

Yreal = data.Yreal;
%Y0 = interp1(data.SvecA1,data.YA1,xi);
Y0 = data.Yreal;
h = mean(Yreal);



%%
%To speed up the ode solving we precompute the intergral I on a square mesh
%We then get specific points of via interplation
% Set up the Intergartion mesh.
    
%Set up mesh
Nxi = N; Neta = 200;
xi_mesh = linspace(-0.5, 0.5-dxi, Nxi);
eta_mesh = linspace(-h, h, Neta+1) ;

deltol = 2*h/(Neta); %distance we offset above the real axis.

% Make mesh
[XI, ETA] = meshgrid(xi_mesh, eta_mesh);
ZETA = XI + 1i*ETA;
% Initialise
Ymat = 0*ZETA;


%%
mytol = 1e-10;


%%

%Need to be cearful when taking paths on real axis in lower half plane
path = [0.05 0.05+0.01*1i -0.05+0.01*1i -0.05... 
    -0.05-0.015*1i 0.05-0.015*1i 0.05-1e-3*1i...
    -0.05-1e-3*1i -0.05-0.015*1i 0.05-0.015*1i 0.05-1e-3*1i... 
    -0.05-1e-3*1i];

 % path = [0.05 0.05+0.01*1i -0.05+0.01*1i -0.05+0.001*1i... 
 %     0+0.001*1i 0-1i*0.01];

%path = [path,path];



omytol = 1e-10;
Ifunc = @(Z) getI(Yreal,Z,N);
myfunc = @(Z,Y) getDY(Y,F,C,Z,Ifunc);
sfun = @(t) eta_mesh(1) + (eta_mesh(end) - eta_mesh(1))*t; %We need sfunc for correc stop conditons 
opts = odeset('Events', @(t,y) steepEvent_time(t, y,sfun,Ifunc,F,C),...
    'RelTol',mytol,...
    'AbsTol', mytol,...
    'Stats','off', ...
    'MaxStep', 6e-4);

Yinit = interp1(xi, Yreal, real(path(1)));

[Svec, out] = get_ode_line(path, Yinit, myfunc, opts);

dX = 1-(Ifunc(Svec)-1i*(sign(imag(Svec))).*out);
dY = myfunc(Svec,out);

J = dX.^2+dY.^2; 

%%
figure(4)
hold on
plot3(real(path),imag(path),min(Yreal)*ones(length(path),1))

plot3(real(Svec),imag(Svec),imag(out),'-')


figure(5)
hold on
plot3(real(Svec),imag(Svec),real(out),'-')
plot3(xi_mesh,0*xi_mesh,Yreal)

hold off

% figure(5)
% hold on
% plot3(real(Svec),imag(Svec),real(J),'o-')
% 

function [value, isterminal, direction] = steepEvent_time(t, y,sfun,Ifunc,F,C)
    %Stop condtion for ode.
    %Stop when the dardient of dy/dx is to steep.
    z = sfun(t);
    % dydt = myfunc(z, y);
    % 
    % value = 1e8   - max(abs(dydt));
    testCon = (F^2/2+C-y)*(1-Ifunc(z));
    value = abs(testCon)-6.8e-3;

    isterminal = 1;             % stop the integration
    direction = -1;             % only trigger when value crosses downward


end
