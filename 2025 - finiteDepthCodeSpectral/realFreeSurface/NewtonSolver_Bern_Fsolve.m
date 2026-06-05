function [Y,C,F,B,M,err,fval]=NewtonSolver_Bern_Fsolve(Y,C,F,B,M,E,iAmp,iEig)

N = length(Y);
L = 0.5;
dxi = 2*L/N;
xi = -L:dxi:L-dxi;
n = [0:N/2-1 0 -N/2+1:-1];
Mk = 1i*n*pi/L;

Yhat = fft(Y);

if iEig==0
    start = [Yhat(1:N/2) F C];%Froude number eigenvalue
    NotEig = B;
elseif iEig==1
    start = [Yhat(1:N/2) B C];%Bond number eigenvalue
    NotEig = F;
end

options = optimoptions('fsolve','Display','iter');
options.MaxIterations = 50;
options.MaxFunctionEvaluations = N*1000;
options.FunctionTolerance=10^(-10);

%Outputing fval which is the neumerical exit stae of the sovler HW 28/07
[sol,ber,fval] = fsolve(@(solution)Steady_Bern_Fsolve(solution,NotEig,E,M,N,n,Mk,dxi,iAmp,iEig),start,options);
err = norm(ber)^2;

YhatHalf = sol(1:N/2);
if iEig==0
    F = sol(N/2+1);%Froude number eigenvalue
elseif iEig==1
    B = sol(N/2+1);%Bond number eigenvalue
end
C = sol(N/2+2);

Yhat = zeros(1,N);
Yhat(1:N/2) = YhatHalf;
Yhat(N/2+1) = 0;
Yhat(N/2+2:end) = fliplr(Yhat(2:N/2));
Y = real(ifft(Yhat));
end


