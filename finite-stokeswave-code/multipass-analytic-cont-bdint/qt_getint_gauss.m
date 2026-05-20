function I = qt_getint_gauss(w, wave)
    % GETINT_GAUSS computes beta integral using a specialized scheme
    %
    %   2015-08: Original implementation of Crew + Trinh
    %   2016-11-30: Finite-depth code for YJL
    r0 = wave.r0;
    
    F = @(phi) get_theta_series(phi, wave.an, wave.r0, wave.N) ...
                .*cot(pi*(phi - w));
    
    I = pvi_gauss(F, w, -0.5, 0.5);
    
    I = reshape(I, size(w));
    
  
                
                
end