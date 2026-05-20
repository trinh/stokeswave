close all


N=50;
s = 0.1;
ten = 2;
grav_coeff = 1;

%for ep = .3
k = 2*pi;
A = 0.1;

%initial guess - last two are mu and B 
x0 = zeros(1,N-1);

for j=1:(N-1)
    x0(j) = ((-1)^(j))*((8*pi*j)/(k))*(A^j);
end
x0(N) = 0.1;
x0(N+1) = 2;


x0(1:N-1) = zeros(1,N-1);
x0(N) = 1;
x0(N+1) = 1;


x0 = est;


options = optimset('Display', 'iter');
fwd = @(x)capsol(x, s, N, ten, grav_coeff);
[est,fval] = fsolve(fwd,x0,options); % Call solver

%%
%integrate up to get x and y. f = delta + ibeta
a=0;
b=1;
M=101;

phi = linspace(a,b,M);

z = phi;
z2 = ones(1,M);


for k=1:N-1
    z = z + ((1i*est(k))/(2*k*pi))*exp(-2*pi*1i*k*phi);
    z2 = z2 + est(k)*(exp(-2*1i*k*pi*phi));
end

x_dash = real(z2);
y_dash = imag(z2);

figure(1);
plot(real(z),imag(z));
axis equal


%{
myeps = double(1 - 1/((abs(z2(1))^2*(abs(z2((M+1)/2))^2))));
disp(['epsilon = ', num2str(myeps)]);
%}
mymu = est(N);
disp(['mu = ', num2str(mymu)]);
