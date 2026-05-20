%{
N=150;
ep = 0.3;
kap = 0;
grav_coeff = 1;
%initial guess - last two are B and mu 
x0 = zeros(1,N-1);
x0(N) = 0.1;
x0(N+1) = 1;
options = optimset('Display', 'iter');
fwd = @(x)myfun_havelock_altcoll(x, ep, N, kap, grav_coeff);
[est,fval] = fsolve(fwd,x0,options); % Call solver
toc
%}

N=150;

syms x
f2 = (1-x)^(-0.5);
test = zeros(1,N);
for k=1:N
    test(k) = (eval(subs(diff(f2,x,(k-1)),0)))/(factorial(k-1));
end

a = test;


%a(1) = 1;
%a(2:N) = est(2:N);

q = 0.33; %reduction factor
tol = 0.5;

n = N;

z = [0, 0.75*1i, 0.85*1i, 0.95*1i];
b = zeros(N,N);
b(1,1:N) = a;

radius = zeros(length(z));
radius(1) = 1/((abs(a(N)))^(1/(N)));

tic
for k=1:(length(z)-1)
    delta = z(k+1)-z(k);
    vector = henriciMatrix(floor(n*q^(k)),floor(n*q^(k-1)), delta)*transpose(b(k,1:floor(n*q^(k-1))));
    b(k+1,1:floor(n*q^(k))) = vector;
    radius(k+1) = (1/(abs(b(k+1,floor(n*q^(k))))))^(1/(floor(n*q^k)));
end
toc


%Now we plot the extensions
figure();
hold on

for k=1:(length(z))
    x = linspace(real(z(k))-tol*radius(k),real(z(k))+tol*radius(k),50);
    y = linspace(imag(z(k))-tol*radius(k),imag(z(k))+tol*radius(k),50);
    [X,Y] = meshgrid(x,y);
    
    f = zeros(50,50);

    for j=1:floor(n*q^(k-1));
        f = f + b(k,j).*(X+1i*Y-z(k)).^(j-1); 
    end

    surf(X,Y,imag(f));
end


[X2,Y2] = meshgrid(linspace(-1,1.5,50),linspace(-1,1.5,50));
f_actual = eval(subs(f2,X2+1i*Y2));
surf(X2,Y2,imag(f_actual));



