%A script to get the analytic contiunation of the leading order surface
%profile.

%This works by solving an ODE for Yxi along verticale strips from the real
%axis up into the upper half plain.

%We requie a data file which contains
%M  'mass'
%F  'Froude number'
%Yreal 'Real solotuion for the free surface'
%Y0     'Solotion of Y on real axis for non primary sheets can be gotten
%callManualPath

clear
close all

%laod data
%data = load("A1realLine.mat");
%data = load("SamDeep.mat");

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

Yreal = data.Yout;
%Y0 = interp1(data.Svec,data.Y1,xi);
Y0 = data.Yout;

h = mean(Yreal);

%%
%Set up mesh
Nxi = N; Neta = 100;
xi_mesh = linspace(-0.5, 0.5-dxi, Nxi);
eta_mesh = linspace(0, 0.2, Neta+1);

% Make mesh
[XI, ETA] = meshgrid(xi_mesh, eta_mesh);
ZETA = XI + 1i*ETA;

%Make sub mesh for calculating
index = 1:1:length(xi_mesh);
Ymat = nan*ZETA(:,index);


%%
disp('Pre calculating I')


%Imat = Copy_of_getI2(Yreal,ZETA,N);
%Ifunc = @(Z) interp2(XI,ETA,Imat,real(Z),imag(Z));

%%
%set up ODE we need to solve to get paths in complex plain.
mytol = 1e-10;
Ifunc = @(Z) getI(Yreal,Z,N);
myfunc = @(Z,Y) getDY(Y,F,C,Z,Ifunc);
sfun = @(t) eta_mesh(1) + (eta_mesh(end) - eta_mesh(1))*t; %We need sfunc for correc stop conditons 
opts = odeset('Events', @(t,y) steepEvent_time(t, y, myfunc,sfun,Ifunc,F,C),...
    'RelTol',mytol,...
    'AbsTol', mytol,...
    'Stats','off', ...
    'MaxStep', 1e-2);
%%
%calculating surface
for count=1:length(xi_mesh)/1
    count
    
    path1 = xi_mesh(index(count))+[eta_mesh(1), eta_mesh(end)]*1i; 
    %path2 = xi_mesh(index(count))+[eta_mesh(Neta/2+1), eta_mesh(1)]*1i;
    Yinit = @(xi0) interp1(xi, Y0, xi0);
    
    [Svec1, y1] =  getOdePath(path1, Yinit(xi_mesh(index(count))),Yreal,F,C,N,M,opts);
    %[Svec2, y2] =  getOdePath(path2, Yinit(xi_mesh(index(count))),Yreal,F,C,N,M);

    %Ymat(:,count) = interp1([flip(imag(Svec2(2:end)))' imag(Svec1)'], [flip(y2(2:end)') y1'], eta_mesh);
    Ymat(:,count) =  interp1(imag(Svec1)',y1',eta_mesh);
end

YmatFull= interp2(XI(:,index),ETA(:,index),Ymat,XI,ETA);


%%
plotSurfaces

%%

function [value, isterminal, direction] = steepEvent_time(t, y, myfunc,sfun,Ifunc,F,C)
    %Stop condtion for ode.
    %Stop when the dardient of dy/dx is to steep.
    z = sfun(t);
    % dydt = myfunc(z, y);
    % 
    % value = 1e8   - max(abs(dydt));
    testCon = (F^2/2+C-y)*(1-Ifunc(z));
    value = abs(testCon) - 1e-3;

    isterminal = 1;             % stop the integration
    direction = -1;             % only trigger when value crosses downward


end
