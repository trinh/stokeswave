function alpha = get_2dcumarc(z)
% GET_CUMARC gets cumulative arclength of a 1D complex curve
%   
%
%   Given a 2D curve via vectors real(z) and imag(z) we find alpha so 
%   that x = x(alpha) and y = y(alpha) where alpha begins at 0 and goes to
%   1, but is scaled according to the segment lengths in between

    W = [real(z(:)) imag(z(:))];    

    %update arclength
    alpha = 0*z;
    alpha(1)=0;
    for j = 2:length(z)
        alpha(j) = norm(W(j,:)-W(j-1,:));
    end

    alpha = cumsum(alpha);    
end

