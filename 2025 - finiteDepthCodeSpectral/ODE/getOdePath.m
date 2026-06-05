function [SvecOut,yOut] = getOdePath(path,Y0,Yreal,F,C,N,M,opts)
%Function to calculate the ode on the fulle space along some path.
h = mean(Yreal);
%% Split path up into sections (Assume striagt up)
if imag(path(2))>2*h
    newPath = [path(1),real(path(1))+i*2*h,path(2)];
    %keyboard
else
    newPath = path;
end
path = newPath;


%% Caes 1 -2h<Imag(z)<2h
%set up ODE we need to solve to get paths in complex plain.
mytol = 1e-10;
Ifunc = @(Z) getI(Yreal,Z,N);
myfunc = @(Z,Y) getDY(Y,F,C,Z,Ifunc);


%%
SvecOut = [path(1)]; yOut = [Y0]; YResPath = [];
for count = 1:length(path)-1

    sfun = @(t) path(count) + (path(count+1) - path(count))*t; %We need sfunc for correc stop conditons 
   

    [Svec, y] = get_ode_line([path(count),path(count+1)], Y0, myfunc,opts);
    SvecOut = [SvecOut; Svec(2:end)];
    yOut = [yOut; y(1:end-1)];
    Y0=y(end);
    
    Ifunc = @(Z) getI(Yreal,Z,N)-2i*interp1(imag(SvecOut),y,imag(Z)-2*h);
    myfunc = @(Z,Y) getDY(Y,F,C,Z,Ifunc);

end







end



