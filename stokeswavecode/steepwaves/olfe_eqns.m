function G = olfe_eqns(b,N)

% b is the vector of coefficients. c = b(N+1). Truncate the series at N

% Shift the elements of b. This is bn' = bn-1
c = b(N+1);
b(2:N+1) = b(1:N);
b(1) = 1;

% Calculate the Bn' (Bn' = Bn-1)
for m=0:N
    sumlimit = N-m;
    B(m+1) = 0;
    for l=0:sumlimit
        B(m+1) = B(m+1) + b(l+1)*b(m+l+1);
    end
end

% Calculate the An'
for n=0:N
    A(n+1) = 0;
    for n2 = 0:n
        A(n+1) = A(n+1) + B(n2+1)*B(n-n2+1);
    end
    A(n+1) = 0.5*A(n+1);
    
    sumlimit2 = N-n;
    for n3 = 1:sumlimit2
        A(n+1) = A(n+1) + B(n3+1)*B(n+n3+1);
    end
end

% Calculate the Fn'
for n=0:N
    F(n+1) = 0;
    for l=0:N
        F(n+1) = F(n+1) + ((6*l+1)*b(l+1))/(9*((2*n+1)^2)-(6*l+1)^2);
    end
    
    F(n+1) = ((18*sqrt(3))/(pi))*F(n+1);
end

% The main equations
for n=0:N
    G(n+1) = (3*n + 2)*A(n+1)- (3*n + 1)*A(n+1) - (F(n+1))/(c^2);
end

end

