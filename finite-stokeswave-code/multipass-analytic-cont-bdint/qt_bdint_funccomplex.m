function [F, theta] = qt_bdint_funccomplex(w, Y, wave, number_of_pass, a)
% QT_BDINT_FUNCCOMPLEX is the nested analytic continuation function
%
%   Y is column vector of size (number_of_pass*length(w))
%       (Standard DE solver will only use sheet_number x 1)
%
%   number_of_pass = 1 => UHP
%   number_of_pass = 2 => UHP, LHP
%   number_of_pass = 3 => crossed axis 2 times 
%
%   a = 1 if the first analytic continuation was into the UHP
%
%   Note the solutions are tau_1 = Y(1) [the desired solution], tau_2 =
%   Y(2) is the first residue switched on, ... tau_n = Y(n) is the last
%   residue (this one is solved independently of any others)
%
%
%   2015-08: Original implementation of Crew + Trinh
%   2016-11-30: Finite-depth code for YJL



if size(Y, 1) ~= number_of_pass
    warning('Y not correct size');    
    % vector inputs...i.e. Y should be (number_of_pass x 1)
    keyboard
end

% Parameters
mu = wave.mu;

theta = NaN*ones(number_of_pass, length(w));
tau = Y;

w = w(:).'; % Convert to row

% Hilbert transform and I-integral are the same size as w
H = 0*w;
I = 0*w;

% For finite depth, we require an additional I integral
Iint = @(phi, f) qt_bdint_Iint(phi, f, wave.an, wave.r0, wave.N);

% Calculate the Hilbert + I-Integral
for k = 1:length(w)
    H(k) = qt_getint_gauss(w(k), wave);
    if wave.r0 == 0
        I(k) = 0;
    else
        I(k) = integral(@(phi)Iint(phi, w(k)), -0.5, 0.5);        
    end
end

% We set the 'bottom' sheet, i.e. does not require recursion
theta(number_of_pass,:) = 1/(1i*(-1)^(number_of_pass+1))*(H + I - tau(number_of_pass,:));

% We have to calculate theta and tau backwards
for j = number_of_pass-1:-1:1    
    summ = 0*w;
    for k = j + 1:number_of_pass
        summ = summ + (-1)^k*2i*theta(k,:);
    end
    theta(j,:) = a/(1i*(-1)^(j+1))*(H + I - tau(j,:) + summ);    
end

% Bernoulli equation is the same
F = -2*pi/mu*sin(theta)./exp(3*tau);
