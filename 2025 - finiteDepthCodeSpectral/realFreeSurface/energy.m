function E=energy(Y,Mk,B,F,n,dxi)
Yhat=fft(Y);
Y1 = real(ifft(Mk.*Yhat));
H = mean(Y);
Mt = 1i*coth(2*pi*H*n);
X1 = 1-real(ifft(Mt.*Mk.*Yhat));
J = X1.^2+Y1.^2;
mass=dxi*sum(Y.*X1);

KE = .5*F^2*dxi*sum(Y.*(X1-1));
CE = B*dxi*sum((J.^(1/2)-X1));
GE = .5*dxi*sum(Y.^2.*X1);

E=KE+CE+GE-0.5*mass^2;


end