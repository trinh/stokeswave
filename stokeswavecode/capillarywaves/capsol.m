function F = capsol(a, s, N, ten, grav_coeff)

%a the vector of unknowns. I.e. a_n coefficients and mu and beta.
%ep given wave "steepness"
%N step size

for j=1:N
    %collocation points
    phi = 0.5*((j-1)/(N-1));
    
    %This is delta + ibeta (without the first term)
    for k=1:N-1
       b(k) = a(k)*exp((-2*k*1i*pi*phi));
    end
    
    %The integral term (gravity)
    for k=1:N-1
       c(k) = (1/(2*k*pi))*(cos(2*k*pi*phi)-1)*a(k);
    end
    
    %The surface tension term
    for k=1:N-1
        d1(k) = -2*pi*k*a(k)*cos(2*pi*k*phi); %beta_phi
        d2(k) = -2*pi*k*a(k)*sin(2*pi*k*phi); %delta_phi
    end
    
    f = 1+sum(b);
    delta = real(f);
    beta = imag(f);
    absval_f = (abs(f))^2;
    
    kap = ten;
    
    tension = ((-kap)/(2*pi*a(N)))*(delta*(sum(d1))-beta*(sum(d2)))/((absval_f)^(3/2));
    grav = 2*pi*(1/a(N))*sum(c);
    energy = 1/(2*absval_f);
    
    F(j) = grav_coeff*grav + energy + tension - a(N+1);
   
end

%{
%steepness condition expressed in terms of ep
for k=1:N-1
    a2(k) = a(k)*exp((-k*1i*pi));
end

F(N+1) = 1 - 1/((abs(1+sum(a(1:N-1)))^2)*(abs(1+sum(a2))^2)) - s;
%}

%in terms of s
for k=1:N-1
    a2(k) = (1/(2*k*pi))*(cos(k*pi)-1)*a(k);
end
F(N+1) = sum(a2) + s;






   
