function I = qt_finite_getint_gauss(w, wave)
    % GETINT_GAUSS computes beta integral using a specialized scheme
    %
    %   2015-08: Original implementation of Crew + Trinh
    %   2016-11-30: Finite-depth code for YJL
    %   2026-05-19: Renaming finite so as to avoid catastrophic problems
    %
    % Note that although this is labeled as "finite" the functional form is
    % the same. It is the integral of:
    %
    %   theta*cot(pi(phi - w)) from w = -0.5 to w = 0.5

   
    
    F = @(phi) getTheta(phi, wave).*cot(pi*(phi - w));
    I = pvi_gauss(F, w, -0.5, 0.5);
    
    I = reshape(I, size(w));

    % Get only the theta index
    function theta = getTheta(phi, wave)        
        [~, theta] = wave.getQTValues(phi);
    end
    
end