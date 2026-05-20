
function [F, theta] = gen_finite_ODE(w, Y, wave, number_of_pass, a, L, UD)

if size(Y, 1) ~= number_of_pass
    warning('Y not correct size');    
    % vector inputs...i.e. Y should be (number_of_pass x 1)
    keyboard
end

% Parameters
mu = wave.mu;

theta = NaN*ones(number_of_pass, length(w));
tau = Y;


w = w(:).'; % Convert paths rows in matrix w

% Hilbert transform and I-integral are the same size as a path in w
H1 = 0*w;
H2 = 0*w;
I = 0*w;


% For finite depth, we require an additional I integral
Iint = @(phi, f) qt_bdint_Iint(phi, f, wave.an, wave.r0, wave.N);

% Calculate the Hilbert + I-Integral
for k = 1:length(w)
    [H1(k),H2(k)] = finite_qt_getint_gauss(w(k), wave);
    if wave.r0 == 0
        I(k) = 0;
    else
        I(k) = integral(@(phi)Iint(phi, w(k)), -0.5, 0.5);        
    end
end

R=H1-H2+I;

%get matrix of coefficients of thetas and taus in each theta eqn
THETACOF=make_theta_matrix(L,UD);
TAUCOF=make_tau_matrix(L,UD);


   
for l=number_of_pass:-1:1

    if L(l)==0
        upperlower = UD(l);
    elseif L(l)>0
        upperlower = -1;
    end
    
        theta(l,:) = upperlower/1i*(R-tau(l,:));
        for i=l:number_of_pass
            theta(l,:) = theta(l,:) + TAUCOF(l,i)*tau(i,:)...
                + THETACOF(l,i)*theta(i,:); 
        end
end
% Bernoulli equation is the same
F = -2*pi/mu*sin(theta)./exp(3*tau);

