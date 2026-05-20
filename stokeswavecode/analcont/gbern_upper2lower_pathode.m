function gp = gbern_upper2lower_pathode(f, g, solseries, phi, gaxis1)
    
    B = solseries.B;
    mu = solseries.mu;    
    
    % To get the value of f in the third quadrant, walk from g(1) upwards
    % and use values from zfree
    % Initial value is gaxis1
    g0 = interp1(phi, gaxis1, real(f));
    % Set the path
    
    
    path = [real(f), real(f) + 1i, -real(f) + 1i, -real(f) - 1i*imag(f)];
    % Get the value
    [ftmp, gtmp, gptmp, zptmp] = gbern_contour(path, solseries, g0);
    p = zptmp(end);
            
    gp = mu./(2*p.*(mu*B + 1i*pi*g)) - p;  
%     disp(['f = ', num2str(f), ', p = ', num2str(p)]);
end