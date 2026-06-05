function out=Steady_Bern_Fsolve(sol,NotEig,E,M,N,n,Mk,dxi,iAmp,iEig)

YhatHalf = sol(1:N/2);
if iEig==0%Froude number eigenvalue
    F = sol(N/2+1);
    B = NotEig;
elseif iEig==1%Bond number eigenvalue
    F = NotEig;
    B = sol(N/2+1);
end
C = sol(N/2+2);

Yhat = zeros(1,N);
Yhat(1:N/2) = YhatHalf;
Yhat(N/2+1) = 0;
Yhat(N/2+2:end) = fliplr(Yhat(2:N/2));

Y=real(ifft(Yhat));
H = mean(Y);
Mt = 1i*coth(2*pi*H*n);
% X = xi-real(ifft(Mt.*Yhat));

Y1 = real(ifft(Mk.*Yhat));
Y2 = real(ifft(Mk.^2.*Yhat));
X1 = 1-real(ifft(Mt.*Mk.*Yhat));
X2 = -real(ifft(Mt.*Mk.^2.*Yhat));

J = X1.^2+Y1.^2;
kap = (X1.*Y2-Y1.*X2)./J.^(3/2);%curvature
ber=(F^2)./(2*J)+Y-B*kap-C-.5*F^2;
ber=real(fft(ber));%Fourier space

mass=dxi*sum(Y.*X1);
if iAmp==1%Amplitude.
    out=[mass-M ber(1:N/2) 100*((max(Y)-min(Y))/2-E)];
end
if iAmp==2%Energy.
    En=energy(Y,Mk,B,F,n,dxi);
    out=[mass-M ber(1:N/2) 1000*(En-E)];
end
end