%Gravity+capillary deep periodic waves

function F = series_bern(a, ep, N, kap, grav_coeff)

%a the vector of unknowns. I.e. a_n coefficients and mu and beta.
%ep given wave "steepness"
%N step size

mu = a(N);
B = a(N+1);


for j=1:N
    %collocation points
    phi = 0.5*((j-1)/(N-1));
    
    %This is delta + ibeta (without the first term)
    dzdf = 1;
    for k=1:N-1
       dzdf = dzdf + a(k)*exp((-2*k*1i*pi*phi));
    end
    
    %The integral term (gravity)
    y = 0;
    for k=1:N-1
       y = y + 1/(2*k*pi)*(cos(2*k*pi*phi)-1)*a(k);
    end
        
    %The surface tension term
%     for i=1:N-1
%         d1(k) = -2*pi*a(k)*cos(2*pi*k*phi); %beta_phi
%         d2(k) = -2*pi*a(k)*sin(2*pi*k*phi); %delta_phi
%     end
    
    
%     delta = real(dzdf);
%     beta = imag(dzdf);
    absval_f = (abs(dzdf))^2;
    
    % tension = ((-kap)/(2*pi*mu))*(delta*(sum(d1)) - beta*(sum(d2)))/((absval_f)^(3/2));
    tension = 0;
    grav = 2*pi/mu*y;
    energy = 1/(2*absval_f);
    
    
    %set the jth equation (Bernoulli condition on free surface)
    F(j) = grav_coeff*grav + energy + tension - B;
end


delta_beta0 = 1 + sum(a(1:N-1));
W0 = 1/(delta_beta0);
delta_beta_end = 1 + sum(a(1:N-1).*exp(-pi*1i*(1:N-1)));
Wend = 1/delta_beta_end;

% F(N+1) = 1 - 1/delta0^2 - ep;
F(N+1) = 1 - abs(W0)^2*abs(Wend)^2 - ep^2;

