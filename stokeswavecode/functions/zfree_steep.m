function [z, zp] = zfree_steep(f, steepsolseries, L)
% ZFREE_steep constructs the series solution using the Havelock steep wave
% method
% M is the number of points for the trapezoid integration
% f should be in the LHP

bn = steepsolseries.bn;
beta = steepsolseries.beta;

[w, ~] = vb_wfree_steep(f, bn, beta);

% Get zp and integrate to get z.
zp = 1./w;

z = 0*f;
for k=1:length(f)
    
    t1 = linspace(0, abs(real(f(k))), L);
    t2 = linspace(0, -imag(f(k)), L);
    
    dgam1 = sign(real(f(k)));
    dgam2 = -1i;
    
    [w_int1, ~] = vb_wfree_steep(sign(real(f(k))).*t1, bn, beta);
    [w_int2, ~] = vb_wfree_steep(-t2.*1i + real(f(k)), bn, beta);
    
    zp_int1 = (1./w_int1)*dgam1;
    zp_int2 = (1./w_int2)*dgam2;
    
    z(k) = trapz(t1, zp_int1) + trapz(t2, zp_int2);    
end

z = transpose(z);
end