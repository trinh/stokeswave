function [Svec, W_mat] = get_ode_line(path, Winit, myfunc,opts)
% GET_ODE_LINE solves an ODE for W(s) along a specified path 

if nargin == 0
    path = exp(1i*linspace(0, 4*pi, 20));
    Winit = 1;
    myfunc = @(t, W) 1./(2*W);
end

Svec = []; % vector for independent variable
W_mat = []; % matrix containing W

%Note the loop only use if we have 2+ points
for count=1:length(path)-1

    % Construct a straight line between each edge in 'points'
    s0 = path(count);
    s1 = path(count+1);
    sfun = @(t) s0 + (s1 - s0)*t;
    ds_dt = s1 - s0;

    % Define an ODE function with alternate parameterisation 
    fwd = @(t, W) ds_dt*myfunc(sfun(t), W);

    [tt, Wsol] = ode23s(fwd,[0 1], Winit, opts);    
    ss = sfun(tt);    
    
    % Extend our solution structures to include newest structure    
    W_mat = [W_mat; Wsol];
    Svec = [Svec; ss];

    %For multi point paths we have to set initial conditions as the last
    %value of W
    Winit = Wsol(end,:);

end

if nargin == 0
    plot3(real(Svec), imag(Svec), real(W_mat));
end


