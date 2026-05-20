function F = olfe_eqns(b,N)

% b is the vector of coefficients. c = b(N+1). Truncate the series at N
% terms. 
% One equation for each unknown

% Shift the elements of b
d(1) = 1;
d(2:N+1) = b;
b=d;
c = b(N+1);

% Calculate all the extra crap

% Create the Bn
for m=1:N+1
    sumlimit = N+1-m;
    B(m) = 0;
    for l=1:sumlimit
        B(m) = B(m) + b(l)*b(m+l);
    end
end

% Create the An
for n=1:N+1
    A(n-1) = 0.5*(sum(B(1:(n-1)).*B(n-1:-1:1))); 
    sumlimit2 = N+1-n;
    for n2 = 1:sumlimit2
        A(n-1) = A(n-1) + B(n2-1)*B(N+1+n2);
    end
end

% Create the Fn
for n=1:N+1
    F(n-1) = 0;
    for n2=1:N+1
        F(n-1) = F(n-1) + (((6*(n2-1)+1)*(b(n2-1)))/(9*((2*(n-1)+1)^2) - (6*(n2-1)+1)^2));
    end
end

% The main equations
for n=1:N+1
    F(n) = (3*(n-1) + 2)*A(n-1)- (3*(n-1) + 1)*A(n) - (F(n-1))/(c^2);
end

end

