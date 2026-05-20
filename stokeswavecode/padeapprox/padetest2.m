% =====================================
%         Load in the solution
F = load('data/series_N101_ep0.7.mat');
solseries = F.solseries;

mu = solseries.mu;
B = solseries.B;
an = solseries.an;
N = solseries.N;
% =====================================


M=100; % Mesh
J=((N-1)/2);
L=J;

x1 = linspace(-5,5,M);
y1 = linspace(-5,5,M);
phi = linspace(-0.5,0.5,M);
psi = linspace(-0.5,0.5,M);
[Phi, Psi] = meshgrid(phi,psi);
[X, Y] = meshgrid(x1,y1);

% Set the coefficients
coefficients(1) = 1;
for k=2:N;
    coefficients(k) = an(k-1);
end

radius = 1/((abs(coefficients(N)))^(1/N));

[P,Q] = padeapprox(coefficients, J, L);

% Pade approximant of z'(f)
for k1=1:M
    for k2=1:M
        numerator = 0;
        denominator = 0;
        for j1=1:J+1
            numerator = numerator + P(J-j1+2)*(exp(-2*pi*1i*(j1-1)*(Phi(k1,k2)+1i*Psi(k1,k2))));
        end
        
        for j2=1:L+1
            denominator = denominator + Q(L-j2+2)*(exp(-2*pi*1i*(j2-1)*(Phi(k1,k2)+1i*Psi(k1,k2))));
        end
       
        zp(k1,k2) = numerator/denominator;
    end
end


% Pade approximant of z(f)

% Set the coefficients of z(f)
coefficients2(1) = 0;
for j=2:N
    coefficients2(j) = ((1i*est(j-1))/(2*pi*(j-1)));
end

[P2,Q2] = padeapprox(coefficients2, J, L);


for k1=1:M
    for k2=1:M
        numerator = 0;
        denominator = 0;
        for j1=1:J+1
            numerator = numerator + P2(J-j1+2)*(exp(-2*pi*1i*(j1-1)*(Phi(k1,k2)+1i*Psi(k1,k2))));
        end
        
        for j2=1:L+1
            denominator = denominator + Q2(L-j2+2)*(exp(-2*pi*1i*(j2-1)*(Phi(k1,k2)+1i*Psi(k1,k2))));
        end
       
        z(k1,k2) = Phi(k1,k2) + 1i*Psi(k1,k2) + numerator/denominator;
    end
end
  


figure(1);
surf(Phi,Psi,imag(z));

figure(2);
surf(Phi,Psi,real(z));




