clear

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
%%
%To speed up the ode solving we precompute the intergral I on a square mesh
%We then get specific points of via interplation
% Set up the Intergartion mesh.
    
%Set up mesh
Nxi = N; Neta = 200;
xi_mesh = linspace(-0.5, 0.5-dxi, Nxi);
eta_mesh = linspace(0, 2*h, Neta+1) ;

deltol = 2*h/(Neta); %distance we offset above the real axis.

% Make mesh
[XI, ETA] = meshgrid(xi_mesh, eta_mesh);
ZETA = XI + 1i*ETA;

%%
Ifunc = @(Z) getI(Yreal,Z,N);

Imat = Ifunc(ZETA);
%%
figure(1)
contourf(XI,ETA,real(Imat))

figure(2)
surf(XI,ETA,real(Imat),'EdgeColor','none')

figure(3)
contourf(XI,ETA,imag(Imat))

figure(4)
surf(XI,ETA,imag(Imat),'EdgeColor','none')

