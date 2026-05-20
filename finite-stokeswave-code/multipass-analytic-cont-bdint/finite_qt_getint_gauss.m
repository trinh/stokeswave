function [I,J] = qt_getint_gauss(w, wave)
    % GETINT_GAUSS computes beta integral using a specialized scheme
    %
    %   2015-08: Original implementation of Crew + Trinh
    %   2016-11-30: Finite-depth code for YJL
    %   2017-01-20: Finite-depth 
    r0 = wave.r0;
    
    F = @(phi) get_theta_series(phi, wave.an, wave.r0, wave.N) ...
                .*cot(pi*(phi - w));
    
            %IF USING THIS TURN I OFF
    FCHECK = @(phi) real(2*(get_tau_series(phi, wave.an, wave.r0, wave.N)...
              -1i*get_theta_series(phi, wave.an, wave.r0, wave.N)).*(exp(-2*pi*...
              1i*phi))./(exp(-2*pi*1i*phi)-exp(-2*pi*1i*w))); 
            
    I = pvi_gauss(F, w, -0.5, 0.5);
    
    I = reshape(I, size(w));
    
    G = @(phi) 2*r0^2*((get_tau_series(phi, wave.an, wave.r0, wave.N) ...
                .*(r0^2-cos(2*pi*(phi - w))))./(r0^4-2*r0^2*cos(2*pi*(phi - w))+1)-...
                (get_theta_series(phi, wave.an, wave.r0, wave.N).*...
                sin(2*pi*(phi - w)))./(r0^4-2*r0^2*cos(2*pi*(phi - w))+1));
      
     GCHECK = @(phi) real( 2*r0^2*(get_tau_series(phi, wave.an, wave.r0, wave.N)...
              +1i*get_theta_series(phi, wave.an, wave.r0, wave.N)).*(exp(-2*pi*...
              1i*phi))./(r0^2*exp(-2*pi*1i*phi)-exp(-2*pi*1i*w)));
     
     
     J = bottom_pvi_gauss(G, w, -0.5, 0.5, r0);
     J = reshape(J, size(w));
                
                
end